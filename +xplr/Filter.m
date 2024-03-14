classdef Filter < xplr.DataOperand
    % function F = filter(header_in[,label])
   
    properties (SetAccess='protected')
        % input: header_in is already a property of the dataOperand mother class
        % operation:
        selection = xplr.SelectionND.empty(1,0);
        % output: header_out is already a property of the dataOperand mother class
    end
    properties (SetObservable, AbortSet)
        slice_fun = @brick.nmean;   % 'nmean', 'mean', 'max', 'min', etc.
        slice_fun_str = 'brick.nmean';
    end
    properties(Dependent, SetAccess='protected', Transient)
        n_sel
        indices
        slice_fun_str_simple
    end
    
    % Setting and updating filter
%     methods (Static)
%         function s = loadobj(obj)
%             s = obj;
%         end
%     end
    methods
        function F = Filter(header_in, label_out)
            % size and header of the input space
            if ~isa(header_in, 'xplr.Header'), error 'first argument must be an xplr.Header object', end
            F.header_in = header_in;
            
            % header of the output space
            if nargin < 2
                if isscalar(header_in)
                    label_out = [header_in.label, ' ROI'];
                else
                    label_out = [brick.strcat({header_in.label}, '-'), ' ROI'];
                end
            end

            % output header is categorical
            if ~isscalar(header_in)
                % header output will be a mere enumeration (no values)
                F.header_out = xplr.Header(label_out, 0);
            elseif header_in.n_column > 0
                % categorical header: we will keep track of values
                F.header_out = xplr.Header(label_out, header_in.sub_labels, cell(0, header_in.n_column));
            elseif header_in.categorical
                % categorical header with no values: keep track of indices
                F.header_out = xplr.Header(label_out, xplr.DimensionLabel('Index', 'numeric'), cell(0, 1));
            else
                % measure header: keep track of values
                F.header_out = xplr.Header(label_out, header_in.sub_labels, cell(0, 1));
            end
            
            % additional information: ROI
            switch F.nd_in
                case 1
                    F.augment_header('ROI1D', 'Selection1D');
                case 2
                    F.augment_header('ROI2D', 'Selection2D');
            end
        end
        function update_selection(F, varargin)
            % function update_selection(F,value)
            % function update_selection(F,'new|all',value[,'label1',header_values1,...])
            % function update_selection(F,'new|all|chg|add',ind,value[,'label1',header_values1,...])
            % function update_selection(F,'chg|add',ind,value)
            % function update_selection(F,'remove|perm',ind)
            % function update_selection(F,'reset')
            
            % input
            if isscalar(varargin)
                if ischar(varargin{1})
                    if strcmp(varargin{1}, 'reset')
                        flag = 'all';
                        value = xplr.SelectionND.empty(1, 0);
                    else
                        error 'only flag ''reset'' can be used without arguments'
                    end
                else
                    flag = 'all';
                    value = varargin{1};
                end
                add_header_info = cell(2, 0);
            else
                flag = varargin{1};
                varargin(1) = [];
                switch flag
                    case {'all', 'new', 'chg', 'add', 'chg&new', 'chg&rm'}
                        if mod(nargin, 2) == 1
                            % ind not specified
                            value = varargin{1};
                            switch flag
                                case 'all'
                                    ind = 1:length(value);
                                case 'new'
                                    ind = F.n_sel + (1:length(value));
                                otherwise
                                    error 'incorrect number of arguments'
                            end
                            add_header_info = reshape(varargin(2:end), 2, []);
                        else
                            % ind specified
                            [ind, value] = deal(varargin{1:2});
                            add_header_info = reshape(varargin(3:end), 2, []);
                        end
                    case {'remove', 'perm'}
                        ind = varargin{1};
                        value = []; % should not be used
                    otherwise
                        error('flag ''%s'' not handled by xplr.filter.update_selection')
                end
            end
            
            % compute indices
            if brick.ismemberstr(flag, {'all', 'new', 'chg', 'add', 'chg&new', 'chg&rm'})
                value = value.compute_indices([F.header_in.n]);
                % filter definition when data is categorical can only be
                % indices based (but not shape based)
                if any([F.header_in.categorical])
                    value = value.convert('indices', [F.header_in.n]);
                end
            end
            
            % update selection
            switch flag
                case 'all'
                    % note that even if value is already equal to
                    % F.selection, we cannot just return, because
                    % add_header_info might bear some changes
                    ind = 1:length(value);
                    F.selection = brick.row(value); % in particular brick.row(...) transform 0x0 array in 1x0 array
                case {'new', 'chg'}
                    F.selection(ind) = value;
                case 'add'
                    if ~isscalar(ind), error 'only one selection at a time can be augmented with flag ''add''', end
                    F.selection(ind) = F.selection(ind).union(value); % automatic union of indices as well
                    flag = 'chg'; % for the notification
                case 'chg&new'
                    F.selection([ind{:}]) = value;
                case 'chg&rm'
                    F.selection(ind{1}) = value;
                    F.selection(ind{2}) = [];
                case 'remove'
                    F.selection(ind) = [];
                case 'perm'
                    F.selection = F.selection(ind);
            end
            
            % update header of output space
            if brick.ismemberstr(flag, {'all', 'new', 'chg', 'chg&new', 'chg&rm'})
                % make a nice name for new 1D selections
                if F.nd_in == 1 && ~isempty(value) && ~any(strcmpi(add_header_info(1, :), 'name'))
                    n = length(value);
                    names = cell(n, 1);
                    for i = 1:n
                        lines = value(i).polygon;
                        sub_names = cell(1, size(lines, 2));
                        for j = 1:length(sub_names)
                            k_start = ceil(lines(1, j));
                            k_stop = floor(lines(2, j));
                            if k_start == k_stop
                                if F.header_in.categorical || F.header_in.is_datetime
                                    try
                                        sub_names{j} = F.header_in.get_item_names{k_start};
                                    catch
                                        sub_names{j} = '(out)';
                                    end
                                else
                                    % measure
                                    [start, scale, unit] = deal(F.header_in.start, F.header_in.scale, F.header_in.unit);
                                    sub_names{j} = sprintf('%.4g%s', start+(k_start-1)*scale, unit);
                                end
                            else
                                if F.header_in.categorical || F.header_in.is_datetime
                                    sub_names{j} = [F.header_in.get_item_names{k_start}, '-', F.header_in.get_item_names{k_stop}];
                                else
                                    % measure
                                    [start, scale, unit] = deal(F.header_in.start, F.header_in.scale, F.header_in.unit);
                                    sub_names{j} = sprintf('%.4g-%.4g%s', start+([k_start k_stop]-1)*scale, unit);
                                end
                            end
                        end 
                        names{i} = brick.strcat(sub_names, ',');
                    end
                    add_header_info(:, end+1) = {'Name', names};
                end
                % ROI
                roi_label = brick.cast(F.nd_in, 'ROI1D', 'ROI2D');
                add_header_info(:, end+1) = {roi_label, value};
                % track header values
                switch flag
                    case 'chg&new'
                        ind_chg_new = [ind{:}];
                    case 'chg&rm'
                        ind_chg_new = ind{1};
                    otherwise
                        ind_chg_new = ind;
                end
                head_value = slice_header(F, ind_chg_new, add_header_info);
                % update header_out
                F.header_out = update_header(F.header_out, flag,ind, head_value);
            else
                F.header_out = update_header(F.header_out, flag, ind);
            end

            % notification
            e = xplr.EventInfo('filter', flag, ind, value);
            notify(F, 'changed_operation', e)
        end
        function set.slice_fun_str(F, fun)
            F.slice_fun_str = fun;
            F.slice_fun = eval(['@', fun]);
        end
        function set.slice_fun(F, fun)
            % input
            if ischar(fun)
                F.slice_fun_str = fun;
            elseif isa(fun, 'function_handle')
                F.slice_fun_str = char(fun);
            else
                error 'incorrect slicing function'
            end
            % some functions are recognized and adapted to be called
            % with the second argument being the dimension
            f_spec = {@std, @nanstd, @brick.nstd, @var, @nanvar, @brick.nvar};
            for i = 1:length(f_spec)
                if isequal(fun, f_spec{i})
                    fun = @(x, d)fun(brick.float(x), 0, d);
                    break
                end
            end
            f_spec = {@min, @max};
            for i = 1:length(f_spec)
                if isequal(fun, f_spec{i})
                    fun = @(x, d)fun(x, [], d);
                    break
                end
            end
            % try it: make sure that when applied to a given dimension the
            % function reduces the test data to an array of same size
            % except 1 in the specified dimension
            test_data = rand(2, 3, 2);
            x = fun(test_data, 2);
            assert(isequal(size(x), [2, 1, 2]))
            % update property
            F.slice_fun = fun;
            % notification
            notify(F, 'changed_operation', xplr.EventInfo('operation'))
        end
        function set_fun(F, fun)
            % convert char to function handle
            if ischar(fun)
                switch fun
                    case {'mean', 'median', 'mode', 'rms'}
                        fun = str2func(fun);
                    case 'max'
                        fun = @(x, dim)max(x, [], dim);
                    case 'min'
                        fun = @(x, dim)min(x, [], dim);
                    case 'std'
                        fun = @(x, dim)std(x, [], dim);
                    case 'var'
                        fun = @(x, dim)var(x, [], dim);
                    otherwise
                        error('unknown function for slicing ''%s''', fun)
                end
            end
            % update slicing function
            F.slice_fun = fun;
            % notification
            notify(F, 'changed_operation', xplr.EventInfo('filter', 'chg', 1:F.n_sel))
        end
        function copy_in(F, obj)
            % do not call updateing methods because there might be
            % additional information in header_out; change manually the
            % needed properties and raise event
            F.slice_fun = obj.slice_fun;
            F.selection = obj.selection;
            F.header_out = obj.header_out;
            e = xplr.EventInfo('filter', 'all', 1:length(F.selection), F.selection);
            notify(F, 'changed_operation', e)
        end
    end
    
    % Get/Set Dependent
    methods
        function n_sel = get.n_sel(F)
            n_sel = length(F.selection);
        end
        function indices = get.indices(F)
            % get data indices from selections
            indices = cell(1, F.n_sel);
            for i = 1:F.n_sel
                indices{i} = F.selection(i).compute_indices([F.header_in.n]).data_ind;
            end
        end
        function str = get.slice_fun_str_simple(F)
            [fun_str, fun_str_simple] = xplr.Filter.menu_fun_options();
            str = F.slice_fun_str;
            idx = brick.find(str, fun_str);
            if ~isempty(idx)
                str = fun_str_simple{idx};
            end
        end
    end
    
    % Slicing
    methods (Access='private')
        function head_value = slice_header(F, ind, add_header_info)
            % columns to which values have been assigned
            n_ind = length(ind);
            head_value = cell(n_ind, F.header_out.n_column);
            ok_column = false(1, F.header_out.n_column);
            
            % tracking of input header values
            if ~isscalar(F.header_in)
                % no tracking when input space has more than one dimensions
            elseif F.header_in.n_column > 0
                % track values: call to a xplr.Header method
                head_value(:, 1:F.header_in.n_column) = F.header_in.track_values(F.indices(ind));
                ok_column(1:F.header_in.n_column) = true;
            elseif F.header_in.categorical
                % categorical with no values: keep track of indices
                head_value(:, 1) = F.indices(ind);
                ok_column(1) = true;
            else
                % measure header: keep track of values
                for i=1:n_ind
                    str = F.header_in.get_item_names(F.indices{ind(i)});
                    if isscalar(str), str = str{1}; else str = brick.row(str); end
                    head_value{i, 1} = str;
                end
                ok_column(1) = true;
            end
            
            % additional values set together with the selections
            if nargin < 3, add_header_info = cell(2, 0); end
            [head_value, affected_columns] = set_add_header_info(F, head_value, add_header_info);
            ok_column(affected_columns) = true;
            
            % put default values in columns which were not set
            for i = find(~ok_column)
                [head_value{:, i}] = deal(F.header_out.sub_labels(i).default_val);
            end
        end
    end
    methods
        function slic = slicing(F, dat, dims, sel_sub_idx)
            % function slic = slicing(F,dat,dims,sel_sub_idx)
            %---
            % slice data in dimension dim according to filter F
            % dat and slic are simple Matlab arrays
            
            % input
            if nargin < 4, sel_sub_idx = 1:F.n_sel; end
            n_sel_slice = length(sel_sub_idx);
            
            % size
            s = size(dat);
            nd_data = max(max(dims), length(s));
            s(end+1:nd_data) = 1;
            
            % initial reshape
            dbef = 1:min(dims) - 1;
            ok = true(1, nd_data);
            ok(dims) = false;
            ok(1:min(dims))=false;
            daft = find(ok); % faster than setdiff
            dat = brick.reshapepermute(dat, {dbef, dims, daft}); % this won't duplicate the array if dims are consecutive
            
            % slicing
            if n_sel_slice == 1
                slic = F.slice_fun(dat(:, F.indices{sel_sub_idx}, :), 2);
            else
                if isa(dat, 'exch.ndSparse')
                    slic = exch.ndSparse.build([], [], [prod(s(dbef)), n_sel_slice, prod(s(daft))]);
                else
                    data_type = brick.switch_case(class(dat), 'double', 'double', 'single');
                    slic = zeros([prod(s(dbef)), n_sel_slice, prod(s(daft))], data_type);
                end
                for i=1:n_sel_slice
                    slic(:, i, :) = F.slice_fun(dat(:, F.indices{sel_sub_idx(i)}, :), 2);
                end
                if isa(slic, 'exch.ndSparse') && (nnz(slic) / numel(slic) > .2)
                    % convert from sparse to full when density becomes more
                    % than 20%
                    slic = full(slic);
                end
            end
            
            % final reshape
            s_final = s;
            s_final(dims(1)) = n_sel_slice;
            s_final(dims(2:end)) = [];
            slic = reshape(slic, [s_final 1]);
        end
    end
    methods (Access='protected')
        function slic = operation_(F, dat, dims)
            % function slic = operation_(F,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            slic = F.slicing(dat, dims);
        end
        function update_operation_(F, x, dims, slice, flag, ind)
            % function update_operation_(F,x,dims,slice,flag,ind)
            %---
            % x and slice are xplr.xdata objects
            
            % slice
            switch flag
                case 'all'
                    slic = F.slicing(x.data, dims);
                case {'new', 'chg'}
                    slic = F.slicing(x.data, dims, ind);
                case {'remove', 'perm'}
                    slic = [];
                case 'chg&new'
                    slic = F.slicing(x.data, dims, [ind{:}]);
                case 'chg&rm'
                    slic = F.slicing(x.data, dims, ind{1});
                otherwise
                    error('flag ''%s'' not handled', flag)
            end
            slice.updateData(flag, dims, ind, slic, F.header_out); % this will trigger automatic notifications
        end
    end
    
    % Link with selection in real world coordinates
    methods
        function selection_world = operation_data_to_space(F)
            header = F.header_in;
            aff = xplr.AffinityND([header.scale], [header.start_num] - [header.scale_num]);
            selection_world = aff.move_selection(F.selection);
        end
        function update_operation_data_to_space(F, WO, e)
            if ~strcmp(e.type, 'filter'), return, end
            [flag, ind, value] = deal(e.flag, e.ind, e.value);
            if ~isempty(value)
                aff = xplr.AffinityND([F.header_in.scale], [F.header_in.start] - [F.header_in.scale]);
                value = aff.move_selection(value);
            end
            switch flag
                case 'all'
                    WO.operation = value;
                case {'new' 'chg'}
                    WO.operation(ind) = value;
                case 'remove'
                    WO.operation(ind) = [];
                case 'perm'
                    WO.operation = WO.operation(ind);
                case 'chg&new'
                    WO.operation([ind{:}]) = value;
                case 'chg&rm'
                    WO.operation(ind{1}) = value;
                    WO.operation(ind{2}) = [];
                otherwise
                    error('flag ''%s'' not handled', flag)
            end
            notify(WO,'changed_operation',xplr.EventInfo('filter', flag, ind, value))
        end
        function update_operation_space_to_data(F, selection_world, e)
            if nargin >= 3
                [flag, ind, value] = deal(e.flag, e.ind, e.value);
            else
                flag = 'all';
                value = selection_world;
                ind = 1:length(value);
            end
            if ~isempty(value)
                aff = xplr.AffinityND(1./[F.header_in.scale], 1 - [F.header_in.start]./[F.header_in.scale]);
                value = value.apply_affinity(aff);
            end
            F.update_selection(flag, ind, value)
        end
    end
    
    % Context menu
    methods (Static)
        function [fun_str, fun_str_simple] = menu_fun_options()
            fun_str = {'brick.nmean', 'brick.nmedian', ...
                'brick.nsum', 'min', 'max', 'std'};
            fun_str_simple = {'mean', 'median', ...
                'sum', 'min', 'max', 'std'};
        end
    end
    methods
        function context_menu(F, m)
            delete(get(m, 'children'))
            [fun_str, fun_str_simple] = xplr.Filter.menu_fun_options();
            brick.propcontrol(F, 'slice_fun_str', ...
                {'menuval', fun_str, [fun_str_simple, {'other...'}]}, ...
                {'parent', m, 'label', 'Filter operation'});
        end
    end
end
