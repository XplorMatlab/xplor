classdef ZoomFilter < xplr.DataOperand
    % function z = zoomfilter(header_in[,zoom[,bin]])
    %---
    % defines zoom_ing, but also binn_ing
   
    properties (SetAccess='protected')
        % filter definition
        zoom = ':'      % ':' or [istart istop] with 1 <= istart <= istop <= header_in.n
        bin = 1
        % output
        indices_in       % data points that will be extracted
        indices_out      % data positions after zoom_ing AND BINING (is equal to indices_in only if bin=1)
    end
    properties (Dependent, SetAccess='protected', Transient)
        zoom_value       % [istart istop]
    end
    
    % Setting and updating filter
    methods
        function z = ZoomFilter(header_in, zoom, bin)
            % invalid zoomfilter
            if nargout == 0, return, end

            % input
            if nargin < 3, zoom = ':'; end
            if nargin < 4, bin = 1; end
            if ~isscalar(header_in)
                z = xplr.ZoomFilter.empty(1, 0);
                for i=1:length(header_in)
                    z(i) = xplr.ZoomFilter(header_in(i));
                end
                return
            end
            
            % input header
            z.header_in = header_in;
            
            % operation definition
            z.zoom = zoom;
            z.bin = bin;
            
            % set indices and output header
            prepare_filter(z, true, true)
        end
        function set_zoom(z, zoom, bin)
            chg_zoom = ~isequal(zoom, z.zoom);
            chg_bin = (nargin >= 3 && bin ~= z.bin);
            if ~chg_zoom && ~chg_bin, return, end
            % check and assign zoom
            if strcmp(zoom, ':')
                z.zoom = zoom;
            else
                if ~isnumeric(zoom) || length(zoom) ~= 2, error 'wrong zoom value', end
                zoom = [max(.5,zoom(1)), min(z.header_in.n+.5,zoom(2))];
                if diff(zoom) <= 0, zoom = ':'; end % invalid zoom -> zoom reset
                z.zoom = zoom;
            end
            % assign bin and update output
            if nargin >= 3
                if bin == 0 || mod(bin,1), error 'binn_ing value must be a positive integer', end
                z.bin = bin;
            end
            prepare_filter(z, chg_zoom, chg_bin) % this will raise 'ChangedOperation' event
        end
        function move_zoom(z, n_step)
            if strcmp(z.zoom, ':'), return, end
            d = diff(z.zoom);
            if n_step > 0
                z2 = min(z.header_in.n + .5, z.zoom(2) + d*n_step);
                z.set_zoom([z2-d, z2])
            else
                z1 = max(.5, z.zoom(1) + d*n_step);
                z.set_zoom([z1, z1+d])
            end
        end
        function set_bin(z, bin)
            if bin == z.bin, return, end
            % check
            if bin == 0 || mod(bin, 1), error 'binn_ing value must be a positive integer', end
            % assign and update output
            z.bin = bin;
            prepare_filter(z, false, true) % this will raise 'ChangedOperation' event
        end
        function copy_in(z, obj)
            z.set_zoom(obj.zoom, obj.bin);
        end
    end
    methods (Access='private')
        function prepare_filter(z, chg_zoom, chg_bin)
            n_in = z.header_in.n;
            [cur_bin, cur_i_out] = deal(z.bin, z.indices_out);
            
            % indices
            if strcmp(z.zoom, ':')
                if z.bin == 1
                    n_out = n_in;
                    z.indices_in = 1:n_in;
                else
                    n_out = floor(n_in/z.bin);
                    z.indices_in = reshape(1:z.bin*n_out, z.bin,n_out);
                end
            else
                if z.bin == 1
                    n_out = fn_coerce(floor(1+diff(z.zoom)), 1, n_in);
                    idx_1 = round(z.zoom(1)+(1+diff(z.zoom)-n_out)/2);
                    z.indices_in = idx_1 + (0:n_out-1);
                else
                    n_out = fn_coerce(floor(1+diff(z.zoom)/z.bin), 1, floor(n_in/z.bin) );
                    idx_1 = round(z.zoom(1) + (1+diff(z.zoom)-n_out*z.bin)/2);
                    idx_1 = fn_coerce(idx_1, 1, n_in-z.bin*n_out + 1);
                    z.indices_in = reshape(idx_1 + (0:z.bin*n_out-1), z.bin,n_out);
                end
            end
            z.indices_out = mean(z.indices_in, 1);
            
            % output header
            head_in= z.header_in;
            if strcmp(z.zoom, ':') && z.bin == 1
                % not any change
                z.header_out = head_in;
            elseif head_in.is_measure
                % if header is a measure, new positions are
                % straightforward to compute
                z.header_out = xplr.Header(head_in.sub_labels, n_out, ...
                    head_in.scale*z.bin, ...
                    head_in.start + (z.indices_out(1)-1)*head_in.scale);
            elseif head_in.n_column == 0
                % no values, keep track of index
                z.header_out = xplr.Header(head_in.label, xplr.DimensionLabel('Index', 'numeric'), num2cell(z.indices_out(:)));
            elseif z.bin == 1
                % no binn_ing: getting values is straightforward
                z.header_out = xplr.Header(head_in.label, head_in.sub_labels, head_in.values(z.indices_in,:));
            else
                % binn_ing
                head_value = z.header_in.trackValues(num2cell(z.indices_in,1));
                z.header_out = xplr.Header(head_in.label, head_in.sub_labels, head_value);
            end
            
            % notifications
            chg_n_out = (n_out ~= length(cur_i_out));
            if chg_zoom
                notify(z, 'ChangedOperation', xplr.EventInfo('zoom',chg_n_out))
            end
            if chg_bin
                notify(z, 'ChangedOperation', xplr.EventInfo('bin'))
            end
            any_chg = chg_n_out || (z.indices_out(1) ~= cur_i_out(1));
            zoom_in = any_chg && (isempty(z.indices_out) ...
                || ((z.bin==1) && ~isempty(cur_i_out) && (z.indices_out(1)>=cur_i_out(1)) && (z.indices_out(end)<=cur_i_out(end)) && (cur_bin==1)));
            if zoom_in
                idx_first = find(cur_i_out==z.indices_out(1), 1, 'first');
                idx_last = find(cur_i_out==z.indices_out(end), 1, 'last');
                idx_rm = [1:idx_first-1, idx_last+1:length(cur_i_out)];
                notify(z, 'ChangedOperation', xplr.EventInfo('filter','remove',idx_rm))
            elseif chg_n_out
                notify(z, 'ChangedOperation', xplr.EventInfo('filter','all'))
            elseif any_chg
                notify(z, 'ChangedOperation', xplr.EventInfo('filter','chg',1:n_out))
            end
        end
    end
    
    % Get Dependent
    methods
        function x = get.zoom_value(z)
            x = z.zoom;
            if strcmp(x, ':'), x = [.5, z.header_in.n+.5]; end
        end
    end
    
    % Slicing
    methods
        function slic = slicing(z, dat, dims, sel_sub_idx)
            % here z can be non-scalar!
            if length(dims)~=length(z), error 'number of dimensions does not match number of points', end
            do_sub_idx = (nargin >= 4);
            
            % size
            s = size(dat);
            nd_data = max(max(dims), length(s));
            s(end+1:nd_data) = 1;
            s_out = s;
            head_out = [z.header_out];
            s_out(dims) = [head_out.n];
            
            % extract sub-data
            subs = substruct('()', repmat({':'}, 1, length(s)));
            no_filt = true;
            for i=1:length(z)
                if ~strcmp(z(i).zoom,':') || z(i).bin ~= 1 || do_sub_idx
                    no_filt = false;
                    ind = z(i).indices_in;
                    if do_sub_idx, ind = ind(:, sel_sub_idx); end
                    subs.subs{dims(i)} = ind; 
                end
            end
            if no_filt
                slic = dat;
            else
                slic = subsref(dat, subs);
            end
            
            % bin
            scur = size(slic);
            scur(end+1:nd_data) = 1; % size before binn_ing
            for i=1:length(z)
                if z(i).bin > 1
                    slic = reshape(slic, [prod(s_out(1:dims(i)-1)), z(i).bin, s_out(dims(i)), prod(scur(dims(i)+1:end))]);
                    slic = nmean(slic, 2);
                end
            end
            slic = reshape(slic, s_out);
        end
    end
    methods (Access='protected')
        function slic = operation_(z, dat, dims)
            % function slic = operation_(z,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            slic = z.slicing(dat, dims);
        end
        function update_operation_(z, x, dims, slice, flag, ind)
            % there is no 'smart' way of updating operation: just do the
            % slicing
            slic = z.slicing(x.data, dims);
            slice.update_data(flag, dims, ind, slic, z.header_out); % this will trigger automatic notifications
        end
    end
    
    % Link with zoom definition in real world coordinates
    methods
        function zoom_world = operation_data_to_space(z)
            if strcmp(z.zoom, ':')
                zoom_world = ':';
            else
                zoom_world = z.header_in.start + (z.zoom-1)*z.header_in.scale;
            end
        end
        function update_operation_data_to_space(z, wo, evnt)
            if ~strcmp(evnt.type,'zoom'), return, end
            wo.operation = z.operation_data_to_space();
            notify(wo, 'ChangedOperation')
        end
        function update_operation_space_to_data(z, world_operation, ~)
            if strcmp(world_operation, ':')
                z.set_zoom(':')
            else
                zoom_idx = 1 + (world_operation - z.header_in.start)/z.header_in.scale;
                z.set_zoom(zoom_idx)
            end
        end
    end
    
    % Tools
    methods
        function idx_1 = orig_to_zoomed(z, idx)
            b = (idx >= z.indices_in(1)) & (idx <= z.indices_in(end));
            idx_1 = idx;
            idx_1(~b) = 0;
            idx_1(b) = 1 + floor((idx(b)-z.indices_in(1))/z.bin);
        end
    end
    
end
