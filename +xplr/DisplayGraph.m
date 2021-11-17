classdef DisplayGraph < xplr.GraphNode
    % This class takes care of conversions between slice and graph coordinates,
	% and of x/y ticks

    % properties
    properties (SetAccess='private') % Meant to be private!
        D % parent xplr.viewdisplay object
        % graphic objects
        xy_ticks 
        separation_lines
        map_lines
        map_face_coloring
    end
    properties (SetObservable, AbortSet=true)
        use_ROI2D_map = true;
        show_separation = false;
        separation_color = [.8, .8, .8];
        show_grid_labels = true;
    end
    % shortcuts
    properties (Dependent, SetAccess='private') % Meant to be private!
        ha
    end
    properties (Dependent, SetAccess='private')
        map_drawing
    end
    % parameters
    properties (SetAccess='private') % Meant to be private!
        x_sep = .2;
        y_sep = .2;
    end
    % pre-computations
    properties (SetAccess='private')
        layout
        zslice_sz    % current size of the zoomed slice
        filling     % how much of the space available for each dimension if filled (vector of values <=1)
    end
    properties (SetAccess='private') % Meant to be private!
        steps
    end
          
    % Constructor, Get dependent
    methods
        function G = DisplayGraph(D)
            % parent xplr.viewdisplay object
            G.D = D;
        end
        function ha = get.ha(G)
            ha = G.D.ha;
        end
        function hp = get.map_drawing(G)
            hp = G.map_lines;
        end
    end
    
    % Pre-computations
    methods (Access='private')
        function [st, sz, fill, x_pair, y_pair] = compute_steps_private(G, org_in)
            % actual steps computation occurs here
            % size and organization
            do_signal = strcmp(G.D.display_mode, 'time courses'); % slight difference for x between time courses and images
            header = G.D.zslice.header;
            sz = [header.n];
            nd = length(sz);
            sz_orig = G.D.slice.sz;
            ok_dim = (sz_orig > 1);
            x_layout = org_in.x;
            y_layout = org_in.y;
            nx = length(x_layout);
            ny = length(y_layout);
            zoom = G.get_zoom('value');              % the 'true' coordinates (i.e. in data before zooming) of space edges
            [idx_offset, bin] = G.get_zoom('off&bin'); % coordinates conversions between original and zoomed data
            
            % zooming
            % It shall be noted that while the zoom specification are
            % arbitrary real values, the data will be cut at integer
            % values, therefore there are two ways to display it:
            % 1 - either shifting it by the mismatch between these real /
            %   integer values (in this case moving the zoom will result in
            %   smooth movements)
            % 2 - either occupying the whole dedicated space without caring
            %   about this mismatch (in this case moving the zoom will
            %   result in step movements)
           
            % Zoom method 1 "continuous":
            % (convert real zoom values into zoomed data coordinates)
            % we have: xtrue   = .5 + offset + (xzoomed-.5)*bin
            % and:     xzoomed = .5 + (xtrue-.5-offset)/bin
            zm1 = .5 + (mean(zoom) - .5 - idx_offset)./bin;   % zoom centers in zoomed data coordinates
            ze1 = diff(zoom)./bin;                       % zoom extents in zoomed data coordinates
            
            % Zoom method 2 "steps":
            zm2 = (sz+1)/2;
            ze2 = sz;
            
            %             % Choose method 2 for internal coordinates anf for "2D grid"
            %             % organization. Choose method 1 otherwise.
            %             [zm ze] = deal(zm1, ze1);
            %             if org_in.xy
            %                 din = [x_layout y_layout];
            %             else
            %                 lastx = find(ok_dim(x_layout),1,'last'); % index of "extern" x dimension
            %                 lasty = find(ok_dim(y_layout),1,'last'); % index of "extern" y dimension
            %                 din = [x_layout(1:lastx-1) y_layout(1:lasty-1)]; % "internal" dimensions
            %             end
            %             zm(din) = zm2(din);
            %             ze(din) = ze2(din);
            
            % Always choose method 2
            [zm, ze] = deal(zm2, ze2);
            
            % (in addition, a small correction is applied for signals,
            % whose edges will be at 1 and n instead of .5 and n+.5 for
            % images and grid cells)
            if do_signal && nx >=1 && sz(x_layout(1)) > 1
                ze(x_layout(1)) = ze(x_layout(1)) - 1;
            end 

            % specific data aspect ratio for some dimensions?
            x_head = header(x_layout);
            y_head = header(y_layout);
            [x_pair, y_pair] = checkpairs(x_head, y_head, do_signal);
            [w, h] = brick.pixelsize(G.D.ha);
            axis_ratio = h/w;

            % G.filling will get values assigned both for 'grid' and
            % 'linear' arrangements
            st = struct;
            [fill, st.x_span, st.y_span] = deal(zeros(1, nd), zeros(1, nx), zeros(1, ny));
            
            % does one dimension have 2D grid organization
            x_avail = 1;
            y_avail = 1; % total available x- and y-span
            if ~isempty(org_in.xy)
                % which dimension
                st.xy_dim = org_in.xy;
                xy_mode = org_in.xy_mode;
                xy_header = header(st.xy_dim);
                
                if strcmp(xy_mode, 'map')
                    % Get the map!
                    idx_roi2d = xy_header.get_column_index('ROI2D');
                    if ~idx_roi2d, error 'no map available for map display mode!'; end
                    rois = [xy_header.values{:,idx_roi2d}]; % xplr.SelectionND object
                    n_roi = length(rois);
                    polys = {rois.polygon};
                    
                    % Coordinates of map sides (use data sizes if available)
                    data_sizes = rois(1).data_sizes;
                    if ~isempty(data_sizes)
                        flip_up_down = true; % ROIs were drawn from an image, put the smaller y-ordinates top (first line in image)
                        range = [.5 data_sizes(1)+.5; .5 data_sizes(2)+.5];
                    else
                        flip_up_down = false; % ROIs were not drawn from an image, keep small y-ordinates bottom
                        % calculate rectangle that contains all ROIs
                        range = zeros(2,2,n_roi);
                        for k = 1:n_roi
                            poly_k = polys{k};
                            range(:,1,k) = min(poly_k,[],2);
                            range(:,2,k) = max(poly_k,[],2);
                        end
                        range = [min(range(:,1,:),[],3) max(range(:,2,:),[],3)];
                        % (add a small gap)
                        range = range + diff(range,1,2)*[-1 1]/20;
                    end
                    
                    % ROI centers
                    centers = zeros(2,n_roi);
                    for k = 1:n_roi
                        centers(:,k) = brick.nmean(polys{k},2);
                    end
                    
                    % Map map sides to the graph sides
                    scale = [1/diff(range(1,:)); 1/diff(range(2,:))];
                    if flip_up_down
                        scale(2) = -scale(2);
                    end
                    offset = - mean(range,2) .* scale;
                    st.xy_offsets = brick.add(offset, brick.mult(scale, centers)); % values between -.5 and .5
                    for k = 1:n_roi
                        polys{k} = brick.add(offset, brick.mult(scale, polys{k}));
                    end
                    st.map_polys = polys;
                    
                    % Space available inside each ROI: choose it such that,
                    % if all regions where non-intersecting, a specific
                    % ratio (<1) of the total space would be covered
                    target_coverage = .25;
                    coverage_per_region = target_coverage / n_roi;
                    st.xy_steps = ones(2,1) * sqrt(coverage_per_region);
                    [x_avail, y_avail] = brick.dealc(st.xy_steps);
                
                    % xy_n_col and xy_n_row will be used for particular
                    % display features such as separations; do as if there
                    % was only one grid cell occupying the full space
                    [st.xy_n_col, st.xy_n_row] = deal(1);

                else
                    % TODO: the lines below are not correct, remove?
                    % % span
                    % [st.x_span(st.xy_dim) st.y_span(st.xy_dim)] = deal(x_avail,y_avail);

                    % determine number of column: what aspect ratio is desired
                    % for the grid elements?
                    n_elem = xy_header.n;
                    if nx && ny && x_pair(end) && y_pair(end)
                        elem_ratio = abs((y_head(end).scale*y_head(end).n)/(x_head(end).scale*x_head(end).n));
                    else
                        elem_ratio = 1; % this value will be only loosely respected
                    end
                    if strcmp(xy_mode,'xy')
                        n_col = brick.coerce(round(sqrt(n_elem*elem_ratio*x_avail/y_avail/axis_ratio)), [1, n_elem]);
                        n_row = ceil(n_elem/n_col);
                    else
                        n_row = brick.coerce(round(sqrt(n_elem/elem_ratio/x_avail*y_avail*axis_ratio)), [1, n_elem]);
                        n_col = ceil(n_elem/n_row);
                    end
                    [st.xy_n_col, st.xy_n_row] = deal(n_col,n_row);

                    % set grid positions
                    st.xy_offsets = zeros(2, n_elem);
                    i0 = (n_col+1)/2;
                    x_step = x_avail/n_col;
                    j0 = (n_row+1)/2;
                    y_step = -y_avail/n_row;
                    st.xy_steps = [x_step, y_step];
                    for k=1:n_elem
                        switch xy_mode
                            case 'xy'
                                [i, j] = ind2sub([n_col, n_row], k);
                            case 'yx'
                                [j, i] = ind2sub([n_row, n_col], k);
                        end
                        st.xy_offsets(:, k) = [(i-i0)*x_step;
                        (j-j0)*y_step];
                    end
                    if sz_orig(st.xy_dim)>1
                        x_avail = x_avail/n_col/(1+G.x_sep);
                        y_avail = y_avail/n_row/(1+G.y_sep);
                    end
                    fill(st.xy_dim) = 1;
                end
            else
                st.xy_dim = [];
                % there will be only a single grid element, which fits the
                % full available space
                st.xy_offsets = zeros(2, 1);
                st.xy_steps = ones(2,1);
                [st.xy_n_col, st.xy_n_row] = deal(1);
            end
            
            % define steps while ensuring aspect ratio for pairs (start
            % from the last dimensions)
            [st.x_offset, st.x_step] = deal(zeros(1, nx)); % offsets will be adjusted so as to keep data point of "center-zoom" coordinates in the middle of the display
            [st.y_offset, st.y_step] = deal(zeros(1, ny));
            fill([x_layout, y_layout]) = 1;
            ix = nx;
            iy = ny;
            while ix > 0 || iy > 0
                % go down to the next pair
                ixnext = find(x_pair(1:ix), 1, 'last');
                iynext = find(y_pair(1:iy), 1, 'last');
                if isempty(ixnext), [ixnext, iynext] = deal(1); end
                % (x)
                for ix=ix:-1:ixnext
                    d = x_layout(ix);
                    st.x_span(ix) = x_avail;
                    st.x_step(ix) = x_avail / ze(d);
                    st.x_offset(ix) = -zm(d)*st.x_step(ix);   % middle of zoom should be placed at the middle of the available space
                    if sz_orig(d)>1, x_avail = st.x_step(ix) / (1+G.x_sep); end  % available x-span for (ix-1)th dimension
                end
                % (y)
                for iy=iy:-1:iynext
                    d = y_layout(iy);
                    st.y_span(iy) = y_avail;
                    st.y_step(iy) = -y_avail / ze(d);        % start from top of the screen (i.e. higher values of y) rather than bottom
                    st.y_offset(iy) = -zm(d)*st.y_step(iy);   % middle of zoom should be placed at the middle of the available space
                    if sz_orig(d) > 1, y_avail = abs(st.y_step(iy)) / (1+G.y_sep); end % available y-span for (iy-1)th dimension
                end
                
                % arrange values to maintain aspect ratio for the pair if
                % there is a pair
                if isempty(ix) || (ix == 1 && ~x_pair(ix)), break, end
                curratio = abs(st.y_step(iy))/st.x_step(ix) * axis_ratio;
                targetratio = abs(y_head(iy).scale/x_head(ix).scale);
                correction = targetratio/curratio;
                if correction > 1
                    % need to reduce x-span
                    d = x_layout(ix);
                    st.x_span(ix) = st.x_span(ix)/correction;
                    st.x_offset(ix) = st.x_offset(ix) + zm(d)*st.x_step(ix)*(1-1/correction);
                    st.x_step(ix) = st.x_step(ix)/correction;
                    fill(d) = 1/correction;
                    x_avail = x_avail/correction;
                elseif correction < 1
                    % need to reduce y-span
                    d = y_layout(iy);
                    st.y_span(iy) = st.y_span(iy)*correction;
                    st.y_offset(iy) = st.y_offset(iy) + zm(d)*st.y_step(iy)*(1-1*correction);
                    st.y_step(iy) = st.y_step(iy)*correction;
                    fill(d) = correction;
                    y_avail = y_avail*correction;
                end
                ix = ix - 1;
                iy = iy - 1;
            end
            
            % available space inside the most interior dimensions (this
            % will be used in two cases: 1) when there are no data
            % dimensions in the x- or y-axis, i.e. available space is 1;
            % 2) for time courses display, as "data value" becomes as a new
            % dimension below the most interior dimension in the y-axis.
            st.x_available = x_avail;
            if do_signal && y_avail == 1
                % signals would occupy the full vertical space, because
                % there are no data dimensions in y or xy location ->
                % get a nicer display by leaving some gaps above and below
                y_avail = 1/(1 + G.y_sep/2);
            end
            st.y_available = y_avail;
        end
    end
    methods
        function layout_id = set_map_mode(G, layout_id)
            % Set map mode to layout if appropriate
            if ~G.use_ROI2D_map, return, end
            xy_dim_id = layout_id.xy;
            if isempty(xy_dim_id), return, end
            xy_header = G.D.zslice.header_by_id(xy_dim_id);
            if isempty(xy_header), return, end
            idx_roi2d = xy_header.get_column_index('ROI2D');
            if ~idx_roi2d, return, end
            % all conditions are satisfied for map mode
            layout_id.xy_mode = 'map';
        end
        function b = map_display(G)
            % beware, map display is not synonymous of map mode! map
            % display is map mode when there are no x and y
            org = G.D.layout;
            b = strcmp(G.D.display_mode, 'image') ...
                && strcmp(org.xy_mode, 'map') && isempty(org.x) && isempty(org.y);
        end
        function any_chg = compute_steps(G)
            % function any_chg = compute_steps(G)
            %---
            % sets G properties x_offset x_step y_offset y_step and tells
            % whether they were changed
            
            % compute steps
            if nargout > 0, prev_steps = G.steps; end
            G.layout = G.D.layout;
            [G.steps, G.zslice_sz, G.filling, x_pair] = compute_steps_private(G, G.layout);
            
            % any change
            if nargout>0
                 any_chg = ~isequal(G.steps, prev_steps);
            end
        end
    end
    
    % Ticks
    methods (Access='private')
        function minimum_spacing = ticks_minimum_spacing(G, k)
            % k = 1 for 'x', 2 for 'y'
            ax_siz = brick.pixelsize(G.ha);
            ax_siz_inch = ax_siz/get(0, 'ScreenPixelsPerInch');
            minimum_spacing_inch = .2; % optimal space between ticks in inches
            
            % target space between ticks
            minimum_spacing = minimum_spacing_inch / ax_siz_inch(k);   % target spacing in axes coordinates
        end
        function step = nice_step(G, min_step)
            % function step = nice_step(G, min_step)
            %---
            % determine nice step for ticks as the smaller "round" number
            % that is more than min_step
            t10 = log10(abs(min_step));
            tests = [1, 2, 5, 10];
            [~, idx] = find(log10(tests) >= mod(t10,1), 1, 'first');
            step = sign(min_step) * 10^floor(t10) * tests(idx);
        end
        function [tick_values, tick_labels] = nice_values_datetime(G, value_start, value_stop, min_sub_step)
            % there will be two sub-steps per step (i.e. one value label
            % every two ticks)
            n_sub_step = 2;
            min_step = min_sub_step * n_sub_step;
            % determine most appopriate display format and steps
            format = xplr.auto_datetime_format(value_start, value_stop, min_step);
            if min_step > days(.5)
                step = days(G.nice_step(days(min_step))); % duration
            elseif min_step > hours(.5)
                step = hours(G.nice_step(hours(min_step)));
            elseif min_step > minutes(.5)
                step = minutes(G.nice_step(minutes(min_step)));
            else
                step = seconds(G.nice_step(seconds(min_step)));
            end
            % work in seconds in any case by substracting day of
            % value_start
            day_start = dateshift(value_start, 'start', 'day');
            value_start = seconds(value_start - day_start); % double
            value_stop = seconds(value_stop - day_start);
            step = seconds(step); % double
            % tick values for all substeps
            sub_step = step / n_sub_step;
            tick_values = sub_step * (ceil(value_start/sub_step) : floor(value_stop/sub_step)); % data coordinates
            do_label = (mod(tick_values, step)==0);
            % convert back to datetime
            tick_values = day_start + seconds(tick_values);
            % tick labels only for steps
            n_tick = length(tick_values);
            tick_labels = cell(1, n_tick);
            t = tick_values(do_label); t.Format = format;
            tick_labels(do_label) = cellstr(char(t));
        end
        function [tick_values, tick_labels] = nice_values(G, value_start, value_stop, min_sub_step)
            % special: datetime
            if isdatetime(value_start)
                [tick_values, tick_labels] = nice_values_datetime(G,value_start, value_stop, min_sub_step);
                return
            elseif isduration(value_start)
                error 'not implemented yet'
            end
            % there will be two sub-steps per step (i.e. one value label
            % every two ticks)
            n_sub_step = 2;
            min_step = min_sub_step * n_sub_step;
            % actual step that will be used: minimal "nice step" that is
            % larger than minstep
            step = G.nice_step(min_step);
            % tick values for all substeps
            sub_step = step / n_sub_step;
            tick_values = sub_step * (ceil(value_start/sub_step) : floor(value_stop/sub_step)); % data coordinates
            % tick labels only for steps
            n_tick = length(tick_values);
            tick_labels = cell(1, n_tick);
            do_label = (mod(tick_values, step)==0);
            % use thousand/million suffixes
            M =  max(abs([value_start value_stop]));
            if M > 1e6 && step >= 1e4
                suffix = 'M';
                mult = 1e-6;
            elseif M > 1e3 && step >= 10
                suffix = 'k';
                mult = 1e-3;
            else
                suffix = '';
                mult = 1;
            end
            if isempty(suffix)
                tick_labels(do_label) = brick.num2str(tick_values(do_label), 'cell');
            else
                tick_labels(do_label) = brick.num2str(tick_values(do_label) *mult, 'cell', ['%g' suffix]);
            end
        end
    end
    methods
        function set.show_grid_labels(G, value)
            G.show_grid_labels = value;
            G.set_ticks()
        end
        function set_ticks(G)
            st = G.steps;
            org = G.layout;
            
            % remove previous xy ticks
            brick.delete_valid(G.xy_ticks)
            G.xy_ticks = [];
            
            % stop if data is too large for being displayed
            if G.D.no_display, return, end
            
            % x and y ticks
            for k = 1:2
                f = brick.cast(k, 'x', 'y');
                d = G.D.active_dim.(f);
                if isempty(d) || ~any(d == org.(f))
                    % note that method setValueTicks, which is normally
                    % called AFTER setTicks, might set yticks later
                    set(G.D.ha, [f 'tick'], [])
                    continue
                end
                
                % header
                head = G.D.zslice.header(d);
                n = head.n;

                % conversion between data coordinates and graph
                do_measure = head.is_measure;
                jf = find(d == G.layout.(f), 1);
                switch f
                    case 'x'
                        f_off = st.xy_offsets(k,1) + st.x_offset(jf) + sum(st.x_offset(jf+1:end) + st.x_step(jf+1:end));
                        f_step = st.x_step(jf);
                    case 'y'
                        f_off = st.xy_offsets(k,1) + st.y_offset(jf) + sum(st.y_offset(jf+1:end) + st.y_step(jf+1:end));
                        f_step = st.y_step(jf);
                end

                % target space between ticks
                minimum_spacing = G.ticks_minimum_spacing(k);
                f_span = brick.cast(k, st.x_span, st.y_span);
                minimum_spacing = minimum_spacing / min(1/f_span(jf), 2); % let this target increase up to a factor of two when dimension occupies only a fraction of the space

                % target space in data coordinates
                minimum_step = minimum_spacing / abs(f_step);
                
                % different display depending on whether header is measure
                % or categorical
                if do_measure
                    [start, scale] = deal(head.start, head.scale);
                    [start, stop] = deal(start, start+(n-1)*scale);
                    [ticks_data, tick_labels] = G.nice_values(start, stop, minimum_step*scale);
                    if isduration(scale)
                        ticks_idx = 1 + seconds(ticks_data-start) / seconds(scale); % data indices coordinates
                    else
                        ticks_idx = 1 + (ticks_data-start) / scale; % data indices coordinates
                    end
                else
                    % ticks for each data point (display only some of
                    % them if there is not enough space for all)
                    tick_labels = brick.row(head.get_item_names());
                    if minimum_step <= 1                     
                        ticks_idx = 1:n;
                    else
                        step = G.nice_step(minimum_step);
                        if strcmp(tick_labels{1}, '1')
                            % it seems that we have a mere enumeration,
                            % use a smart step
                            ticks_idx = step:step:n;
                        else
                            % make both the first and last appear
                            ticks_idx = 1:step:n;
                            if ticks_idx(end) ~= n
                                ticks_idx(end) = n;
                            end
                        end
                        tick_labels = tick_labels(ticks_idx);
                    end
                end
                tick = f_off + ticks_idx*f_step;
                
                % set ticks!
                if strcmp(f, 'x')
                    set(G.ha, 'xtick', tick, 'xticklabel', tick_labels)
                else
                    set(G.ha, 'ytick', fliplr(tick), 'yticklabel', fliplr(tick_labels))
                end
            end
            
            % grid
            if ~isempty(st.xy_dim) && G.show_grid_labels
                d = st.xy_dim;
                
                % header
                head = G.D.zslice.header(d);
                n = head.n;

                % display labels for each xy grid cell
                tick_labels = brick.row(head.get_item_names());
                if isempty(st.y_span)
                    row_height = st.y_available;
                else
                    row_height = st.y_span(end);
                end
                xy = brick.add([0; row_height/2], st.xy_offsets);
                rotation = brick.switch_case(st.xy_n_col<10, 0, 20);
                G.xy_ticks = gobjects(1, n);
                for i=1:n
                    G.xy_ticks(i) = text(xy(1,i), xy(2,i), tick_labels{i}, ...
                        'parent', G.ha, 'hittest', 'off', ...
                        'Interpreter', 'none', ...
                        'horizontalalignment', 'center', 'verticalalignment', 'baseline', ...
                        'rotation', rotation);
                end
                if d == G.D.color_dim
                    G.color_grid_ticks()
                end
            end
            
            % show separations
            G.draw_separations()
            
        end
        function color_grid_ticks(G)
            st = G.steps;
            if isempty(st.xy_dim), return, end
            d = st.xy_dim;
            head = G.D.zslice.header(d);
            if d == G.D.color_dim
                colors = head.get_color();
                for i=1:head.n
                    set(G.xy_ticks(i), 'color', colors(i, :))
                end
            else
                set(G.xy_ticks, 'color', 'default')
            end
        end
        function set_value_ticks(G)
            % do we show value ticks?
            if ~strcmp(G.D.display_mode, 'time courses') ...
                    || (~isempty(G.D.active_dim.y) && ~isequal(G.D.active_dim.y, G.steps.xy_dim)) ...
                    || G.D.no_display
                ylabel(G.D.ha, '')
                return
            end
            
            % clip values available?
            if isempty(G.D.grid_clip)
                error 'programming: set_value_ticks should be called after viewdisplay.updateDisplay, so grid_clip should be set'
            end
            sz = strict_size(G.D.grid_clip, 1+G.D.nd);
            sz(1) = [];
            
            % enough space on y-axis to show values?
            org = G.D.layout;
            st = G.steps;
            minimum_spacing = G.ticks_minimum_spacing(2);
            if ~isempty([org.y st.xy_dim]), minimum_spacing = minimum_spacing/2; end
%             if st.y_available < target_spacing
%                 set(G.D.ha,'ytick',[])
%                 ylabel(G.D.ha,'')
%                 return
%             end
            
            % do clipping ranges and baselines differ along 'horizontal'
            % dimensions (i.e. dimensions at location x, mergeddata or
            % xy)?
            h_aligned_dims = [org.merged_data, org.x(2:end)];
            h_aligned_dims(G.D.slice.sz(h_aligned_dims)==1) = [];
            y_dims = org.y;
            if isempty(st.xy_dim)
                xy_n_show = 1;
                xy_on_single_column = false;
            elseif strcmp(org.xy_mode, 'map')
                % show value ticks only for the first map element
                xy_n_show = 1;
                xy_on_single_column = false;
                h_aligned_dims(end+1) = G.steps.xy_dim;
            else
                % show one set of value ticks per row of the xy grid
                xy_n_show = st.xy_n_row;
                xy_on_single_column = (st.xy_n_col==1);
                if xy_on_single_column
                    y_dims(end+1) = G.steps.xy_dim;
                elseif ~isempty(st.xy_dim)
                    h_aligned_dims(end+1) = G.steps.xy_dim;
                end
            end
            same_clip_row = ~any(ismember(h_aligned_dims, G.D.clipping.independent_dim));
            if ~same_clip_row
                % not same clipping for different cells aligned
                % horizontally: not possible to display value ticks
                set(G.D.ha, 'ytick', [])
                ylabel(G.D.ha, '')
                return
            end
            same_baseline = isempty(G.D.signals_baseline) || isempty (h_aligned_dims);

            % remove these dimensions and get back to original signals'
            % clipping range by adding baseline if possible
            grid_clip = subsref_dim(G.D.grid_clip, 1+h_aligned_dims, ones(1, length(h_aligned_dims)));
            if ~isempty(G.D.signals_baseline) && same_baseline
                signals_baseline = subsref_dim(G.D.signals_baseline, h_aligned_dims, ones(1, length(h_aligned_dims)));
                grid_clip = brick.add(grid_clip, shiftdim(signals_baseline,-1));
            end

            % now grid_clip is nonsingleon only in the y_dims dimensions;
            % permute dimensions (put these dimensions first)
            grid_clip = permute(grid_clip, [1, 1+y_dims, 1+setdiff(1:G.D.nd, y_dims)]);
            if xy_on_single_column
                grid_clip = reshape(grid_clip, [2, prod(sz(org.y)), xy_n_show]);
            else
                grid_clip = reshape(grid_clip, [2, prod(sz(org.y))]);
                grid_clip = repmat(grid_clip, [1, 1, xy_n_show]);
            end
            
            % build y_tick and y_tick_values
            ny = prod(sz(org.y));
            [y_tick, y_tick_values, y_tick_labels] = deal(cell([ny, xy_n_show]));
            for k_row = 1:xy_n_show
                for k_y = 1:ny
                    % vertical center of the row
                    y_idx = brick.row(brick.indices(sz(org.y), k_y, 'g2i'));
                    switch org.xy_mode
                        case {'xy' 'map'}
                            y_row_offset = st.xy_offsets(2, 1+(k_row-1)*st.xy_n_col);
                        case 'yx'
                            y_row_offset = st.xy_offsets(2, k_row);
                        case ''
                            y_row_offset = 0;
                    end
                    yoffset = y_row_offset + sum(st.y_offset) + sum(st.y_step .* y_idx);
                    % tick values
                    clip_k = grid_clip(:, k_y, k_row);
                    if any(isnan(clip_k)), continue, end % happens when data itself consists only of NaNs
                    clip_extent = diff(clip_k);
                    minimum_sub_step = minimum_spacing * clip_extent/st.y_available;
                    [tick_values, tick_labels] = G.nice_values(clip_k(1), clip_k(2), minimum_sub_step);
                    y_scale = st.y_available / clip_extent;
                    clip_center = mean(clip_k);
                    y_tick_values{k_y, k_row} = tick_values;
                    y_tick{k_y, k_row} = yoffset + (tick_values-clip_center) * y_scale;
                    y_tick_labels{k_y, k_row} = tick_labels;
                end
            end
            % y_offset are descending, so read y_tick in reverse order to
            % have only increasing values
            y_tick = [y_tick{end:-1:1}];
            y_tick_values = [y_tick_values{end:-1:1}];
            y_tick_label = [y_tick_labels{end:-1:1}];
            % if baselines were not the same (so clip values could not be
            % corrected), indicate that clip values are relative to
            % baselines
            if ~same_baseline
                [y_tick_label{y_tick_values==0}] = deal(['(' G.D.clipping.align_signals ')']);
                idx_p = (y_tick_values>0 & ~brick.isemptyc(y_tick_label));
                y_tick_label(idx_p) = brick.map(@(str)['+' str], y_tick_label(idx_p), 'cell');
            end

            % set ticks
            set(G.D.ha, 'ytick', y_tick, 'yticklabel', y_tick_label)
            ylabel(G.D.ha, 'values')
                
            
        end
    end

    % Separations marks for 'external' dimensions and map display
    methods
        function set.use_ROI2D_map(G, value)
            G.use_ROI2D_map = value;
            G.D.reset_display()
        end
        function set.show_separation(G, value)
            G.show_separation = value;
            if G.map_display()
                % patch is created any way; just show or not the edges
                G.update_separation_color()
            else
                G.draw_separations()
            end
        end
        function set.separation_color(G, value)
            G.separation_color = value;
            update_separation_color(G)
        end
        function update_separation_color(G)
            if G.show_separation
                color = G.separation_color;
            else
                color = 'None';
            end
            if G.use_ROI2D_map
                set(G.map_lines, 'EdgeColor', color)
            else
                set(G.separation_lines, 'Color', color)
            end
        end
        function draw_separations(G)
            % Draw either the separations between all grid elements, or the
            % regions when using map display
            
            % no 'smart update', always delete all existing separation
            % lines and redisplay new ones
            brick.delete_valid(G.separation_lines)
            brick.delete_valid(G.map_lines)
            [G.separation_lines, G.map_lines] = deal([]); 
            
            % no display?
            no_display = ~G.show_separation && ~G.map_display();
            if no_display, return, end
            
            % some properties
            org = G.D.layout;
            nd = G.D.zslice.nd;
            sz = G.D.zslice.sz;
            st = G.steps;
            
            % Map or Grid?
            if ~isempty(org.xy) && strcmp(org.xy_mode, 'map')
                % Map
                do_map_display = G.map_display();
                
                % vertices: easy
                vertices = [st.map_polys{:}]';

                % faces: we need to cut ROI into polygon subregions
                % separated by NaNs (that's a bit difficult)
                n_roi = length(st.map_polys);
                n_face_per_roi = zeros(1, n_roi);
                faces = cell(1, n_roi);
                if do_map_display
                    face_coloring = cell(1, n_roi);
                end
                idx_offset = 0;
                n_points = brick.itemlengths(st.map_polys);
                n_max = max(n_points);
                for k = 1:n_roi
                    poly_k = st.map_polys{k};
                    idx_nan = find(any(isnan(poly_k),1));
                    idx_sub = [0 idx_nan; idx_nan n_points(k)+1];
                    n_p_sub = diff(idx_sub) - 1;
                    idx_sub(:, n_p_sub<3) = []; % we want only faces with at least 3 vertices
                    n_sub = size(idx_sub, 2);
                    faces_k = NaN(n_sub, n_max);
                    for i = 1:n_sub
                        faces_k(i, 1:n_p_sub(i)) = idx_sub(1,i)+1:idx_sub(2,i)-1;
                    end
                    faces{k} = idx_offset + faces_k;
                    idx_offset = idx_offset + n_points(k);
                    n_points(k) = max(n_p_sub);
                    if do_map_display
                        face_coloring{k} = zeros(n_sub, n_roi);
                        face_coloring{k}(:, k) = 1;
                    end
                end
                faces = cat(1, faces{:});
                faces = faces(:, 1:max(n_points));
                if do_map_display
                    G.map_face_coloring = cat(1, face_coloring{:});
                end
                
                % display
                if G.show_separation
                    edge_color = G.separation_color;
                else
                    edge_color = 'None';
                end
                lines = patch('Faces', faces, 'Vertices', vertices, ...
                    'FaceColor', 'none', 'EdgeColor', edge_color, ...
                    'Parent', G.D.ha, 'HitTest', 'off');
                G.map_lines = lines;
            else
                % Grid
                lines = {};
                
                % x
                for k = 2:length(org.x)
                    d = org.x(k);                      % dimension for which we will draw vertical lines
                    d_out = [org.x(k+1:end), st.xy_dim];  % other dimensions more external than di in x location
                    n = sz(d);
                    if n == 1, continue, end
                    sz_out = ones(1, nd);
                    sz_out(d_out) = sz(d_out);
                    sz_out(nd+1) = st.xy_n_col;
                    n_out = prod(sz_out);
                    x_pos = zeros(n-1, n_out);
                    for k_out = 1:n_out
                        % indices of external dimensions
                        ijk = repmat(brick.indices(sz_out, k_out, 'g2i'), [1, n-1]);
                        % indices for dimension d
                        ijk(d, :) = 1.5:n-.5;
                        % x-positions of n-1 lines
                        x_pos(:, k_out) = ...
                            sum(brick.add(st.x_offset(k:end)', brick.mult(ijk(org.x(k:end), :), st.x_step(k:end)')), 1) ...
                            + st.xy_offsets(1) + (ijk(nd+1, :)-1) * st.xy_steps(1);
                    end
                    level = ~isempty(st.xy_dim) + length(org.x) - k + 1;
                    color = 1 - (1-G.separation_color) / 2^level;
                    lines{end+1} = brick.lines('x', x_pos(:), G.D.ha, ...
                        'color', color, 'hittest', 'off'); %#ok<AGROW>
                end

                % y
                for k = (1+strcmp(G.D.display_mode, 'image')) : length(org.y)
                    d = org.y(k);                      % dimension for which we will draw vertical lines
                    d_out = [org.y(k+1:end), st.xy_dim];  % other dimensions more external than di in x location
                    n = sz(d);
                    if n == 1, continue, end
                    sz_out = ones(1,nd);
                    sz_out(d_out) = sz(d_out);
                    sz_out(nd+1) = st.xy_n_row;
                    n_out = prod(sz_out);
                    y_pos = zeros(n-1, n_out);
                    for k_out = 1:n_out
                        % indices of external dimensions
                        ijk = repmat(brick.indices(sz_out, k_out, 'g2i'), [1, n-1]);
                        % indices for dimension d
                        ijk(d, :) = 1.5:n-.5;
                        % y-positions of n-1 lines
                        y_pos(:, k_out) = ...
                            sum(brick.add(st.y_offset(k:end)', brick.mult(ijk(org.y(k:end), :), st.y_step(k:end)')), 1) ...
                            + st.xy_offsets(2) + (ijk(nd+1, :)-1) * st.xy_steps(2);
                    end
                    level = ~isempty(st.xy_dim) + length(org.x) - k + 1;
                    color = 1 - (1-G.separation_color) / 2^level;
                    lines{end+1} = brick.lines('y', y_pos(:), G.D.ha, ...
                        'color', color, 'hittest', 'off'); %#ok<AGROW>
                end

                % xy
                if ~isempty(st.xy_dim)
                    x_pos = -.5 + (1:st.xy_n_col-1) / st.xy_n_col;
                    y_pos = -.5 + (1:st.xy_n_row-1) / st.xy_n_row;
                    lines = [lines brick.lines(x_pos, y_pos, G.D.ha, ...
                        'color', G.separation_color)];
                end
                
                % single vector of graphic handles
                lines = brick.map(@brick.row, lines, 'array'); 
                G.separation_lines = lines;
            end

            % put lines below other graphic elements
            uistack(lines, 'bottom')
            
        end
    end
    
    % Coordinates conversions
    methods (Access='private')
        function [sub_dim, ijk0, mode, invertible] = conversion_options(G, np, varargin)
            % Options (name/value pairs) for slice/graph conversions:
            % - 'invertible'  [default false] if set to true, exterior
            %               coordinates are rounded, this make the
            %               conversion invertible by calling slice_to_graph
            % - 'sub_dim'    [default all dims] dimensions in ijk for which
            %               we perform the conversion; other dimensions
            %               will be assigned to the fixed default values in
            %               ijk0
            % - 'ijk0'      [required if sub_dim is set] default values for
            %               dimensions where no conversion is requested
            % - 'mode'      value is 'point' or 'vector' [default]

            p = inputParser;
            p.addParameter('sub_dim', [], @isnumeric)
            p.addParameter('ijk0', [], @isnumeric)
            p.addParameter('mode', 'point', @(s)ismember(s, {'point', 'vector'}))
            p.addParameter('invertible', false, @islogical)
            parse(p, varargin{:})
            s = p.Results;
            [sub_dim, ijk0, mode, invertible] = ...
                deal(s.sub_dim, s.ijk0, s.mode, s.invertible);
            if isempty(sub_dim)
                sub_dim = 1:G.D.nd;
            elseif ~isempty(ijk0) && size(ijk0, 2) == 1 && np>1
                ijk0 = repmat(ijk0, [1, np]);
            end

        end
    end
    methods
        function [zoom, bin] = get_zoom(G, varargin)
            % function zoom = get_zoom(G[,dim][,'value|effective|indices'])
            % function [offset bin] = get_zoom(G[,dim],'off&bin')
            %---
            % Get zoom value for specified dimensions. 
            % Different modes affect the returned value:
            % - 'value' [default]   returns the zooming value specified in
            %               the zoom filter
            % - 'effective' zooming value, after taking into account the
            %               extra space due to not completely filling the
            %               available space (to preserve data aspect ratio)
            % - 'indices'   minimal and maximal displayed data indices
            % - 'displaylimit'   minimal and maximal slice coordinates
            %               being displayed (same as 'indices' for time
            %               courses display, otherwise extends by +/- .5)
            % - 'off&bin'   returns offset (i.e. number of data points from
            %               the beginning which are not shown) and binning
            %               value
            
            % input
            dim = 1:G.D.nd;
            mode = 'value';
            for i=1:length(varargin)
                a = varargin{i};
                if isnumeric(a)
                    dim = a;
                else
                    mode = a;
                end
            end
            
            % output
            z_filters = G.D.zoom_filters(dim);
            switch mode
                case {'value', 'effective'}
                    zoom = cat(1, z_filters.zoom_value)';
                    if strcmp(mode,'effective')
                        zoom = brick.add(mean(zoom), brick.mult([-.5; .5], diff(zoom)./G.filling(dim)));
                    end
                case {'indices', 'displaylimit'}
                    zoom = zeros(2, length(dim));
                    for i=1:length(dim), zoom(:, i) = z_filters(i).indices_out([1, end]); end
                    if strcmp(mode, 'displaylimit') 
                        if strcmp(G.D.display_mode, 'time courses') && ~isempty(G.layout.x)
                            external_x_dim = (dim ~= G.layout.x(1));
                            zoom(:, external_x_dim) = brick.add(zoom(:, external_x_dim), [-.5; .5]);
                        else
                            zoom = brick.add(zoom, [-.5; .5]);
                        end
                    end
                case 'off&bin'
                    offset = zeros(1, length(dim));
                    for i=1:length(dim), offset(i) = z_filters(i).indices_in(1) - 1; end
                    zoom = offset;          % first output: offset
                    bin = [z_filters.bin];   % second output: bin
            end
        end
        function xy = zslice_to_graph(G, ijk, varargin)
            % function xy = zslice_to_graph(G,ijk[,options...])
            %---
            % Input:
            % - ijk         index coordinates in the zslice data
            % - options (name/value pairs): see xplr.DisplayGraph.conversion_options
            %
            % Output:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            %
            % See also xplr.DisplayGraph.conversion_options

            org = G.layout;
            st = G.steps;
            
            % Input points
            np = size(ijk, 2);
            [sub_dim, ijk0, mode, invertible] = conversion_options(G, np, varargin{:});
            if strcmp(mode, 'vector')
                error 'case not handled yet'
            elseif invertible
                warning 'zslice_to_graph conversion is always invertible, no need to use ''invertible'' flag!'
            end
            if length(sub_dim) < G.D.nd
                % do not consider all dimensions
                % - for ignored dimensions that are "more exterior" than
                %   dimensions in sub_dim, we must choose some fixed value
                %   (we choose 1 by default)
                % - ignored dimensions that are "more interior" than 
                %   dimensions in sub_dim can be totally ignored
                if isempty(ijk0), ijk0 = ones(G.D.nd, np); end
                ijk_ = ijk0;
                if size(ijk, 1) == length(sub_dim)
                    ijk_(sub_dim, :) = ijk;
                elseif size(ijk, 1) == G.D.nd
                    ijk_(sub_dim, :) = ijk(sub_dim, :);
                else
                    error('expected %i or %i number of dimensions, but entry points have %i', length(sub_dim), G.D.nd, size(ijk,1))
                end
                ijk = ijk_;
                org_x_ok = find(ismember(org.x, sub_dim), 1, 'first'):length(org.x); % dimensions on x layout to consider; can be empty
                org_y_ok = find(ismember(org.y, sub_dim), 1, 'first'):length(org.y); % dimensions on x layout to consider; can be empty
            else
                org_x_ok = 1:length(org.x);
                org_y_ok = 1:length(org.y);
            end
            
            % "exterior" dimensions must be rounded
            do_round = true(1, G.D.nd);
            if ~isempty(org_x_ok), do_round(org.x(org_x_ok(1))) = false; end
            if ~isempty(org_y_ok), do_round(org.y(org_y_ok(1))) = false; end
            ijk(do_round, :) = round(ijk(do_round, :));
            
            x = sum(brick.add(st.x_offset(org_x_ok)', brick.mult(ijk(org.x(org_x_ok), :), st.x_step(org_x_ok)')), 1);
            y = sum(brick.add(st.y_offset(org_y_ok)', brick.mult(ijk(org.y(org_y_ok), :), st.y_step(org_y_ok)')), 1);
            xy = [x; y];
            if ~isempty(st.xy_dim)
                xy_idx = ijk(st.xy_dim, :);
                inside = (xy_idx > 0 & xy_idx <= size(st.xy_offsets, 2));
                xy(:, inside) = xy(:, inside) + st.xy_offsets(:, xy_idx(inside));
                % points outside of graph
                xy(:, ~inside) = NaN;
            end
        end
        function ijk = graph_to_zslice(G, xy, varargin)
            % function ijk = graph_to_zslice(G,xy,options...)
            %---
            % Input:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            % - options (name/value pairs): see xplr.DisplayGraph.conversion_options
            %
            % Output:
            % - ijk         index coordinates in the zslice data
            %
            % See also xplr.DisplayGraph.conversion_options

            % Input points
            if isvector(xy), xy = xy(:); end
            np = size(xy, 2);
            ijk = ones(G.D.nd, np);
            
            % Parse options
            [sub_dim, ijk0, mode, invertible] = conversion_options(G, np, varargin{:});
            if length(sub_dim) < G.D.nd && isempty(ijk0)
                % Define ijk0 using the first point
                ijk0 = graph_to_zslice(G, xy(:, 1), 'invertible', true);
                if np == 1, ijk = ijk0(sub_dim); return, end
            end
            
            % If mode is 'vector', we cannot operate in xy dims, and
            % operate at most on one x and one y dims
            org = G.layout;
            if strcmp(mode, 'vector')
                ok = ~any(ismember(sub_dim, org.xy)) ...
                    && sum(ismember(sub_dim, org.x)) <= 1 ...
                    && sum(ismember(sub_dim, org.y)) <= 1;
                if ~ok
                    error 'vector conversion not possible in graph_to_zslice for this set of dimensions'
                end
            end          
                        
            % xy
            st = G.steps;
            sz = G.zslice_sz;
            if strcmp(mode, 'point') && ~isempty(st.xy_dim)
                % take advantage on the fact that the grid spans the full
                % axis
                d = st.xy_dim;
                x = brick.coerce(.5 + xy(1, :), .001, .999); % 0 = left edge, 1 = right edge
                y = brick.coerce(.5 - xy(2, :), .001, .999); % 0 = top edge,  1 = bottom edge
                icol = .5 + x * st.xy_n_col;
                irow = .5 + y * st.xy_n_row;
                if ismember(d, sub_dim)
                    if d == G.layout.xy
                        ijk(d, :) = icol + st.xy_n_col * round(irow-1);
                    else
                        ijk(d, :) = irow + st.xy_n_row * round(icol-1);
                    end
                    ijk(d, :) = min(ijk(d, :), sz(d) + .4999);
                else
                    ijk(d, :) = ijk0(d, :);
                end
                xy = xy - st.xy_offsets(:, round(ijk(d, :)));
            end
            
            % x 
            x = xy(1, :);
            x_layout = org.x;
            for ix = length(x_layout):-1:1
                d = x_layout(ix);
                if strcmp(mode, 'point')
                    x = x - st.x_offset(ix);
                end
                if ismember(d, sub_dim)
                    ijk(d, :) = x / st.x_step(ix);
                else
                    ijk(d, :) = ijk0(d, :);
                end
                if strcmp(mode, 'point')
                    x = x - round(ijk(d, :)) * st.x_step(ix);
                end
            end
            
            % y
            y = xy(2, :);
            y_layout = org.y;
            for iy = length(y_layout):-1:1
                d = y_layout(iy);
                if strcmp(mode, 'point')
                    y = y - st.y_offset(iy);
                end
                if ismember(d, sub_dim)
                    ijk(d, :) = y / st.y_step(iy);
                else
                    ijk(d, :) = ijk0(d, :);
                end
                if strcmp(mode, 'point')
                    y = y - round(ijk(d, :)) * st.y_step(iy);
                end
            end
            
            % we want an output that can be invertible by calling
            % zslice_to_graph, this means that we should not give the
            % conversion "per dimension" but in a global fashion where
            % dimensions "exterior" to the dimensions of interest are
            % rounded
            if invertible
                do_round = true(1, length(ijk));
                if length(sub_dim) < G.D.nd
                    do_round(org.x(find(ismember(org.x, sub_dim), 1, 'first'))) = false; % do not round the first dimension of interest on the x location
                    do_round(org.y(find(ismember(org.y, sub_dim), 1, 'first'))) = false; % do not round the first dimension of interest on the x location
                else
                    if ~isempty(x_layout), do_round(x_layout(1)) = false; end
                    if ~isempty(y_layout), do_round(y_layout(1)) = false; end
                end
                ijk(do_round) = round(ijk(do_round));
            end
            
            % not all dimensions?
            if length(sub_dim) < G.D.nd
                ijk = ijk(sub_dim, :);
            end
        end
        function zijk = slice_to_zslice(G, ijk, do_vector, sub_dim)
            if nargin<3, do_vector = false; end
            [idx_offset, bin] = G.get_zoom('off&bin');
            if nargin>=4
                idx_offset = idx_offset(sub_dim);
                bin = bin(sub_dim);
            end
            if do_vector
                zijk = brick.div(ijk, bin(:));
            else
                zijk = brick.div(brick.subtract(ijk, idx_offset(:)) - .5, bin(:)) + .5;
            end
        end
        function ijk = zslice_to_slice(G, zijk, do_vector, sub_dim)
            if nargin<3, do_vector = false; end
            [idx_offset, bin] = G.get_zoom('off&bin');
            if nargin>=4
                idx_offset = idx_offset(sub_dim);
                bin = bin(sub_dim);
            end
            if do_vector
                ijk = brick.mult(zijk, bin');
            else
                ijk = brick.add(idx_offset', .5 + brick.mult(zijk-.5, bin'));
            end
        end
        function xy = slice_to_graph(G, ijk, varargin)
            % function xy = slice2graph(G,ijk[,options...])
            %---
            % Input:
            % - ijk         index coordinates in the zslice data
            % - options (name/value pairs): see xplr.DisplayGraph.conversion_options
            %
            % Output:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            %
            % See also xplr.DisplayGraph.codnversionOptions

            % Input points
            np = size(ijk, 2);
            [sub_dim, ijk0, mode, invertible] = conversion_options(G, np, varargin{:});
            if strcmp(mode, 'vector')
                error 'case not handled yet'
            elseif invertible
                warning 'zslice_to_graph conversion is always invertible, no need to use ''invertible'' flag!'
            end
            if length(sub_dim) < G.D.nd
                if isempty(ijk0)
                    ijk_ = zeros(G.D.nd, np);
                else
                    ijk_ = ijk0;
                end
                if size(ijk, 1) == length(sub_dim)
                    ijk_(sub_dim, :) = ijk;
                elseif size(ijk, 1) == G.D.nd
                    ijk_(sub_dim, :) = ijk(sub_dim, :);
                else
                    error('expected %i or %i number of dimensions, but entry points have %i', length(sub_dim), G.D.nd, size(ijk,1))
                end
                ijk = ijk_;
            end
            
            % first convert from slice to zoomed slice
            zijk = G.slice_to_zslice(ijk, strcmp(mode, 'vector'));
            
            % then convert to graph coordinates
            if length(sub_dim) < G.D.nd && isempty(ijk0)
                % let zslice_to_graph estimate zijk0
                xy = zslice_to_graph(G, zijk, 'mode', mode, 'sub_dim', sub_dim);
            else
                % ijk and zijk values in other dimensions than sub_dim have
                % already be assigned using ijk0
                xy = zslice_to_graph(G, zijk, 'mode', mode);
            end
        end
        function ijk = graph_to_slice(G, xy, varargin)
            % function ijk = graph_to_slice(G,xy,options...)
            %---
            % Input:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            % - options (name/value pairs): see xplr.DisplayGraph.conversion_options
            %
            % Output:
            % - ijk         index coordinates in the slice data
            %
            % See also xplr.DisplayGraph.conversion_options
            
            % if ijk0 argument, convert it to zoomed slice coordinates
            np = size(xy, 2);
            [sub_dim, ijk0, mode, invertible] = conversion_options(G, np, varargin{:});
            if ~isempty(ijk0)
                zijk0 = G.slice_to_zslice(ijk0);
            else
                zijk0 = [];
            end

            % coordinates in zoomed slice
            zijk = G.graph_to_zslice(xy, 'mode', mode, ...
                'sub_dim', sub_dim, 'ijk0', zijk0, 'invertible', invertible);
            
            % convert to before zooming
            ijk = G.zslice_to_slice(zijk, strcmp(mode,'vector'), sub_dim);
        end
	end

	% Specialized position functions
	methods
        function M = get_transform(G, ijk)
            % function M = get_transform(G,ijk,ybase)
            %---
            % Matrix transformation to place curve/image at ijk data
            % coordinates; note that only coordinates not belonging to the
            % curve/image will be taken into account.
            % 
            % For 'image' display_mode, M will transform indices
            % 1:size(im,1) and 1:size(im,2) to accurate pixel positions.
            % For 'time courses' display_mode, M will transform indices
            % 1:length(x) to accurate x-ordinates, and data values 0 and 1
            % respectively to the bottom and top y-ordinates of the
            % available space.
            % 
            % Input:
            % - ijk     nd * npoint array
            
            st = G.steps;
            org = G.layout;

            % Initialize matrix
            n_transform = size(ijk, 2);
            M = repmat(eye(4), [1, 1, n_transform]);
            
            % Scale & offset
            % (x)
            if isempty(org.x)
                % no data dimension on x-axis: only 1 data point for time
                % courses or image display, which must be positionned in
                % the center of the available space
                x_scale = st.x_available;
                % index 1 must be positionned at x-ordinate 0, i.e. x_offset + 1*x_scale = 0
                x_offsets = -x_scale * ones(1, n_transform);
            else
                x_scale = st.x_step(1);
                x_offsets = brick.add(sum(st.x_offset), sum(brick.mult(brick.column(st.x_step(2:end)), ijk(org.x(2:end), :)), 1));
            end
            % (y)
            switch G.D.display_mode
                case 'image'
                    % not possible to have negative scaling in the
                    % hgtransform matrix
                    % -> orienting the images downward will be achieved by
                    % inverting y coordinates at the stage of patch
                    % creation, i.e. will be -1:-1:-size(im,2)
                    if isempty(org.y)
                        % no data dimension on y-axis: similar to above
                        y_scale = st.y_available;
                        % index -1 must be positionned at y-ordinate 0, i.e. y_offset + 1*y_scale = 0
                        y_offsets = y_scale * ones(1, n_transform);
                    else
                        y_scale = abs(st.y_step(1));
                        y_offsets = brick.add(sum(st.y_offset), sum(brick.mult(brick.column(st.y_step(2:end)), ijk(org.y(2:end), :)), 1));
                    end
                case 'time courses'
                    % in this case, the transformation will apply not on
                    % data indices in some dimension, but on data values
                    % -> orient these values upward, and transform them
                    % such that [0 1] fills the available space
                    y_scale = st.y_available;
                    y_offsets = brick.add(sum(st.y_offset), sum(brick.mult(brick.column(st.y_step(1:end)), ijk(org.y(1:end), :)), 1) );
                    % value .5 must be positionned at y-ordinates y_offsets, i.e. y_offsets + .5*y_scale = 0
                    y_offsets = y_offsets - y_scale/2; 
                otherwise
                    error 'invalid display mode'
            end
            % (add offsets for the xy grid)
            xy_offset = [x_offsets; y_offsets];
            if ~isempty(st.xy_dim)
                xy_offset = xy_offset + st.xy_offsets(:, ijk(st.xy_dim, :)); 
            end
            % (set transform matrix)
            M(1, 1, :) = x_scale;
            M(2, 2, :) = y_scale; 
            M(1:2, 4, :) = xy_offset;
        end
        function [bottom_left, siz] = sub_axes_position(G, ijk)
            % function [bottom_left, siz] = sub_axes_position(G,ijk)
            %---
            % Input:
            % - ijk     nd * npoint array
            % 
            % Output:
            % - pos     4 * npoint array
            
            np = size(ijk, 2);
            st = G.steps;

            % Position of sub-axes centers
            % (x)
%             if isempty(org.x)
%                 % no data dimension on x-axis: only 1 data point for time
%                 % courses or image display, which must be positionned in
%                 % the center of the available space
%                 x_offsets = zeros(1,np); 
%             else
                x_offsets = brick.add(sum(st.x_offset(2:end)), sum(brick.mult(brick.column(st.x_step(2:end)), ijk(org.x(2:end), :)), 1) );
%             end
            % (y)
            switch G.D.display_mode
                case 'image'
                    y_offsets = brick.add(sum(st.y_offset(2:end)), sum(brick.mult(brick.column(st.y_step(2:end)), ijk(org.y(2:end), :)), 1));
                case 'time courses'
                    y_offsets = brick.add(sum(st.y_offset(1:end)), sum(brick.mult(brick.column(st.y_step(1:end)), ijk(org.y(1:end), :)), 1));
                otherwise
                    error 'invalid display mode'
            end
            % (add offsets for the xy grid)
            xy_offsets = [x_offsets; y_offsets];
            if ~isempty(st.xy_dim)
                xy_offsets = xy_offsets + st.xy_offsets(:, ijk(st.xy_dim, :)); 
            end
            
            % We are done!
            siz = repmat([st.x_span(1); st.y_span(1)], [1, np]);
            bottom_left = xy_offsets - siz/2;
        end
        function pos = label_position(G, dim, org_in)
            % function pos = label_position(G,d[,org_in])
            
            dim = G.D.slice.dimension_number(dim);
        
            % steps
            if nargin == 3
                st = compute_steps_private(G, org_in);
            else
                org_in = G.layout;
                st = G.steps;
            end
            
            % label positions
            n = length(dim);
            pos = zeros(1, n);
            for i=1:n
                d = dim(i);
                if ~isempty(org_in.x) && d == org_in.x(end)
                    pos(i) = st.xy_offsets(1, 1);
                elseif any(d == org_in.x)
                    ix = find(d == org_in.x, 1);
                    pos(i) = st.xy_offsets(1, 1) + sum(st.x_offset(ix+1:end) + st.x_step(ix+1:end));
                    if pos(i)<-.5, pos(i) = pos(i) + sum(st.x_step(ix+1:end)); end % first grid element is more than half-outside
                elseif ~isempty(org_in.y) && d == org_in.y(end)
                    pos(i) = st.xy_offsets(2, 1);
                elseif any(d == org_in.y)
                    iy = find(d == org_in.y, 1);
                    pos(i) = st.xy_offsets(2, 1) + sum(st.y_offset(iy+1:end) + st.y_step(iy+1:end));
                    if pos(i) > .5, pos(i) = pos(i) + sum(st.y_step(iy+1:end)); end % first grid element is more than half-outside
                elseif d == st.xy_dim
                    pos(i) = 0;
                else
                    pos(i) = 0;
                end
            end
        end
		function [polygon, center] = selection_mark(G, dim, sel)
            % Create the polygon to display corresponding to a given
            % selection. This is a complex function as it handles many
            % different cases whether the selection is 1D or 2D, which
            % dimensions the selection applies to, and where they are
            % located.
            
            org = G.layout;
            
			% checks
			nd = length(dim);
            dim = G.D.slice.dimension_number(dim);
			if sel.nd ~= nd, error 'selection has incorrect number of dimensions', end
            
            % default polygon is empty (no display)
            polygon = nan(2, 1); 
            center = nan(2, 1); % out of display

			switch nd
				case 1
					lines = sel.polygon; % 2*n array: set of lines
					n_line = size(lines, 2);
                    % remove lines that are completely out of current view
					zoom = G.get_zoom(dim, 'value');
                    lines(:, lines(1, :) > zoom(2) | lines(2, :) < zoom(1)) = [];
                    if isempty(lines), return, end
                    % lines spanning beyond the left or right side
                    beyond_left = lines(1, :) < zoom(1);
                    beyond_right = lines(2, :) > zoom(2);
                    % clip lines to current view
                    lines(1, beyond_left) = zoom(1);
                    lines(2, beyond_right) = zoom(2);
                    % convert from slice to zslice coordinates
                    lines = G.slice_to_zslice(lines, false, dim);
                    % display selections as rectangles (for 'x' and 'y'
                    % locations), or as more complex polygon (for 'xy')
                    st = G.steps;
                    dim_location = G.D.layout_id.dim_locations{dim};
					if ismember(dim_location, {'x', 'y'})
						% convert from zslice to graph coordinates:
						% ignore dimensions that are more internal than dim
						% take value 1 for dimensiont that are more external than dim
						dim_layout = org.(dim_location);
						idx_dim = find(dim_layout == dim, 1);
						switch dim_location
							case 'x'
								lines = sum(st.x_offset(idx_dim:end)) + lines*st.x_step(idx_dim) + sum(st.x_step(idx_dim+1:end));
							case 'y'
								lines = sum(st.y_offset(idx_dim:end)) + lines*st.y_step(idx_dim) + sum(st.y_step(idx_dim+1:end));
						end
                        if ~isempty(st.xy_dim)
                            graph_dim = brick.switch_case(dim_location, 'x', 1, 'y', 2);
                            lines = lines + st.xy_offsets(graph_dim, 1); 
                        end
                        % construct polygon as union of rectangles
                        polygon = cell(1, 2*n_line - 1);
                        for i = 1:n_line
                            switch 2*beyond_left(i) + beyond_right(i)
                                case 0
                                    % segment within view: full rectangle
                                    polygon{2*i - 1} = [lines([1, 2, 2, 1, 1], i)'; -.5, -.5, .5, .5, -.5];
                                case 1
                                    % 'rectangle' open on the right side
                                    polygon{2*i - 1} = [lines([2, 1, 1, 2], i)'; -.5, -.5, .5, .5];
                                case 2
                                    % 'rectangle' open on the left side
                                    polygon{2*i - 1} = [lines([1, 2, 2, 1], i)'; -.5, -.5, .5, .5];
                                case 3
                                    % 'rectangle' open on both sides: 2
                                    % lines
                                    polygon{2*i - 1} = [lines([1, 2], :)', NaN, lines([2, 1], i)'; ...
                                        -.5, -.5, NaN, .5, .5];
                            end
                        end
                        [polygon{2:2:end}] = deal([NaN; NaN]);
                        polygon = [polygon{:}];
                        center = [mean(lines(:)); 0];
                        
                        % invert coordinates if dim location is 'y'
                        if strcmp(dim_location, 'y')
                            polygon = polygon([2, 1], :);
                            center = center([2, 1]);
                        end
                    elseif strcmp(dim_location, 'xy')
                        % fancy display of selections in grid !
                        hh = abs(st.xy_steps(2)) *.46; % half height of the frame around a single grid element
                        polygon = cell(1, 2*n_line - 1);
                        for i = 1:n_line
                            % corners: convert from zslice to graph coordinates
                            r_line = round(lines(:, i) + [1; -1]*.01);
                            c = st.xy_offsets(:, r_line); % 2*2, i.e. x/y * start/stop
                            c(1, :) = c(1, :) + st.xy_steps(1) * (lines(:, i) - r_line)';
                            % sub-polygon
                            single_row = (diff(c(2, :)) == 0);
                            if single_row
                                % the easy case: a simple rectangle
                                % spanning a single line in the grid
                                polygon{2*i - 1} = [c(1,[1, 2, 2, 1, 1]); c(2,[1, 1, 2, 2, 1]) + [-1, -1, 1, 1, -1]*hh];
                            elseif diff(r_line) < st.xyn_col
                                % two non-intersecting rectangles on two successive lines
                                polygon{2*i-1} = ...
                                    [c(1,1), .5*ones(1,2), c(1,[1 1]), NaN, -.5, c(1,[2, 2]), -.5*ones(1,2); ...
                                    (c(2,[1, 1])+hh), (c(2,[1, 1])-hh), (c(2,1)+hh), NaN, (c(2,[2, 2])+hh), (c(2,[2, 2])-hh), (c(2,2)+hh)];
                            else
                                % more difficult: spanning multiple lines
                                hhc = abs(st.xy_steps(2)) - hh;
                                polygon{2*i-1} = ...
                                    [c(1,1), .5*ones(1,3), c(1,2)*ones(1,2), -.5*ones(1,3), c(1,1)*ones(1,2); ...
                                    (c(2,[1, 1])+hh), (c(2,1)-hhc), (c(2,[2, 2])+hhc), (c(2,[2, 2])-hh), (c(2,2)+hhc), (c(2,[1, 1])-hhc), c(2,1)+hh];
                            end
                        end
                        [polygon{2:2:end}] = deal([NaN; NaN]);
                        polygon = [polygon{:}];
                        
                        % center: will be better positioned if we average
                        % after conversion from indices to graph positions
                        center = mean(G.slice_to_graph(sel.data_ind, 'sub_dim', dim), 2);
                    else
                        error 'not implemented yet'
                        center = [brick.nmean(polygon(1, :)), brick.nmean(polygon(2, :))];
                    end
                case 2
                    % somehow simpler because only the x,y configuration is
                    % allowed
                    % get polygon in slice coordinates
                    poly_slice = sel.polygon;
                    center_slice = [mean(poly_slice(1, 1:end-1)); mean(poly_slice(2, 1:end-1))]; % remove the last point as it repeats the first one
                    % restrict to the part that is visible within the
                    % current zoom
                    poly_slice = G.visible_polygon(poly_slice, dim);
                    display_limit = G.get_zoom(dim, 'displaylimit');
                    if any(center_slice < display_limit(1, :)' | center_slice > display_limit(2, :)')
                        center_slice(:) = NaN;
                    end
                    % convert to graph
                    polygon = G.slice_to_graph(poly_slice, 'sub_dim', dim);
                    center = G.slice_to_graph(center_slice, 'sub_dim', dim);
				otherwise
					error 'case not handled yet'
			end

        end
        function output = visible_polygon(G, poly_slice, dim)
                % function output = visible_polygon(G, poly_slice, dim);
                %---
                % @param poly_slice: nd*np list of points
                % @param dim: 1*nd list of dimensions in which the
                % polygon is defined
                % @return output: nd*np' list of points: part of the
                % polygon that is within the zoom limits
                % - points that are outside the zoom limits will have their
                % values replaced by NaNs
                % - intermediary points will be inserted exactly at the
                % limit

                % sizes
                [nd, np] = size(poly_slice);
                assert(length(dim) == nd)

                % get zoom limits in these dimensions
                zoom_slice_values = G.get_zoom(dim, 'displaylimit');

                % set the ouput to zeros (they will be set to one if one of the
                % dimension if it's out of display)
                polygon_is_out_of_display = false(1, np);
                for dimension = 1:nd  
                   % is equal to one if is out of limits of the zoom or if the
                   % previous value was already 1
                    polygon_is_out_of_display = poly_slice(dimension, :) < zoom_slice_values(1, dimension) | poly_slice(dimension, :) > zoom_slice_values(2, dimension) | polygon_is_out_of_display;
                end

                % "Boundary points" will be inserted where there are
                % connections between a displayed point and a point out of
                % display.
                % Example:
                % polygon_is_out_of_display = [0  1  1  0  1  0  0  1]
                % boundary_dir =            [1  0 -1  1 -1  0  1]
                % boundary_prev =            1,    3, 4, 5,    7
                % boundary_next =            2,    4, 5, 6,    8
                % insertion_indexes =       [2     5  7  8    10] 
                boundary_dir = diff(polygon_is_out_of_display);
                boundaries_indexes = find(boundary_dir ~= 0);
                n_boundaries = size(boundaries_indexes, 2);
                boundaries_values = nan(nd, n_boundaries);
                insertion_indexes = zeros(1, n_boundaries);

                % calculate intermediate points between points displayed
                % and points not displayed
                number_of_points_added = 0;
                for k = 1:n_boundaries
                    boundary_prev = boundaries_indexes(k);
                    boundary_next = boundary_prev + 1;

                    % which of the two points is inside / outside
                    if(boundary_dir(boundary_prev) == 1)
                        [inside_point, outside_point] = deal(boundary_prev, boundary_next);
                    else
                        [inside_point, outside_point] = deal(boundary_next, boundary_prev);
                    end

                    % scan dimensions to determine what portion of the
                    % initial segment should be hidden.
                    biggest_ratio = 0;
                    for dimension = 1:nd
                        outside_to_inside = poly_slice(dimension, inside_point) - poly_slice(dimension, outside_point);
                        outside_to_limit_min = zoom_slice_values(1, dimension) - poly_slice(dimension, outside_point);
                        outside_to_limit_max = zoom_slice_values(2, dimension) - poly_slice(dimension, outside_point);

                        % portion of initial segment must be hidden only if
                        % both values have the same sign (if they have
                        % opposite sign, this means that 'outside_point' is
                        % inside the zoom limit for this dimension)
                        V = [outside_to_limit_max, outside_to_limit_min];
                        if ~any(diff(sign(V)))
                            biggest_ratio = max(biggest_ratio, min(abs(outside_to_limit_min), abs(outside_to_limit_max))/abs(outside_to_inside));
                        end
                    end

                    % use this ratio to define boundary point
                    boundaries_values(:, k) = poly_slice(:, outside_point) + (poly_slice(:, inside_point) - poly_slice(:, outside_point))*biggest_ratio;

                    % insertion position
                    insertion_indexes(k) =  boundary_next + number_of_points_added;
                    number_of_points_added = number_of_points_added + 1;
                end

                % replace values by NaN for points that must not be
                % displayed
                poly_slice(:, polygon_is_out_of_display) = NaN;

                % insert boundary points
                output = ones(nd, np + size(boundaries_values, 2));
                output(:, setdiff(1:end, insertion_indexes)) = poly_slice;
                output(:, insertion_indexes) = boundaries_values;
        end
        function sel_slice = selection_to_slice(G, dim, sel_ax)
            % More or less the symmetric of the previous function: convert
            % selection definition in graph coordinates to slice
            % coordinates in the dimensions dim
            % selectionnd object -> use appropriate method for affinity.
            % Works currently only with 2D selections.
            
            if length(dim)~=2, error 'number of dimensions must be 2', end
            
            % use the first point as the origin, work in the zslice to
            % avoid difficulties due to binning
            xy0 = sel_ax.shapes(1).points(:, 1);
            zijk0 = round(G.graph_to_zslice(xy0));

            % infer the affinity matrix graph->zslice
            % (first the linear part)
            xy_test = brick.add(xy0, [0, 1, 0; 0, 0, 1]);
            zijk_test = G.graph_to_zslice(xy_test, 'sub_dim', dim, 'ijk0', zijk0);
            linear_part = [zijk_test(:, 2) - zijk_test(:, 1), zijk_test(:, 3) - zijk_test(:, 1)];
            % (then the offset)
            offset = zijk_test(:, 1) - linear_part * xy0;
            
            % now the affinity matrix graph->slice
            [idx_offset, bin] = G.get_zoom('off&bin');
            linear_part = diag(bin(dim)) * linear_part;
            offset = idx_offset(dim)' + .5 + (offset - .5) .* bin(dim)';
            
            % construct affinitynd object
            affinity = xplr.AffinityND(linear_part,offset);
            
            % use the selectionnd method
            sel_slice = sel_ax.apply_affinity(affinity, G.D.slice.sz(dim));
        end
        function [zoom_dim, zslice_zoom, clip_zijk, clip_zoom] = zoom_in_rect(G, rect)
            % determine in which dimension to zoom and zoom values
            %
            % Input:
            % - rect    nd*2 array, first and second point selected by the
            %           mouse
            %
            % Ouput:
            % - zoom_dim    dimension(s) in which to zoom in
            % - zoom        zoom values in this(these) dimension(s)
            
            org = G.layout;
            st = G.steps;
            sz = G.zslice_sz;
            
            % by default, no clip zoom
            [clip_zijk, clip_zoom] = deal([]);
            
            % start with the most external dimensions to determine whether
            % they are valid for zooming in; if not, go to more internal
            % dimensions
            
            % xy
            if ~isempty(st.xy_dim)
                % correct the rectangle: flip y to go from top to bottom,
                % coerce to within the graph, go from top-left to
                % bottom-right corner
                rect_mod = rect; 
                rect_mod(2, :) = -rect_mod(2, :);
                rect_mod = brick.coerce(rect_mod, -.5, .5);
                rect_mod = [min(rect_mod, [], 2) max(rect_mod, [], 2)];
                
                % coordinates of rectangle in grid
                grid_size = [st.xy_n_col; st.xy_n_row];
                rect_grid_pos = round(.5 + brick.mult(rect_mod + .5, grid_size));
                rect_in_cell = rect_mod + .5 - brick.div(rect_grid_pos - .5, grid_size);
                col = rect_grid_pos(1, :);
                row = rect_grid_pos(2, :);
                d = st.xy_dim;
                if d == org.xy
                    idx = (row-1) * st.xy_n_col + col;
                else
                    idx = (col-1) * st.xy_n_row + row;
                end
                idx = min(idx, sz(d));
                
                % is there at least one grid cell that fits completely
                % within the rectangle?
                cell_size = zeros(2, 1);
                if isempty(st.x_span)
                    cell_size(1) = st.x_available;
                else
                    cell_size(1) = st.x_span(end);
                end
                if isempty(st.y_span)
                    cell_size(2) = st.y_available;
                else
                    cell_size(2) = st.y_span(end);
                end
                rect_grid_inside = rect_grid_pos + ...
                    [(rect_in_cell(:, 1) > -cell_size/2), -(rect_in_cell(:, 2) < cell_size/2)];
                if d == org.xy
                    idx_inside = (rect_grid_inside(2, :)-1) * st.xy_n_col + rect_grid_inside(1, :);
                else
                    idx_inside = (rect_grid_inside(1, :)-1) * st.xy_n_row + rect_grid_inside(2, :);
                end
                idx_inside = brick.coerce(idx_inside, 1, sz(d));
                if diff(idx_inside) >= 0
                    % yes we are covering at least one cell, so we have a
                    % zoom in the "grid" dimension
                    zoom_dim = d;
                    zslice_zoom = idx_inside + [-.4999 .4999];
                    return
                end
                
                % no valid zoom in the "grid" dimension; we will consider
                % zoom in more internal dimensions only if the rectangle
                % covers, even partially, no more than one grid element
                rect_grid_covered = rect_grid_pos + ...
                    [(rect_in_cell(:, 1) > cell_size/2), -(rect_in_cell(:, 2) < -cell_size/2)];
                if any(diff(rect_grid_covered, 1, 2))
                    return
                end
                
                % get rectangle coordinates inside the selected grid
                % element
                col = rect_grid_covered(1, 1);
                row = rect_grid_covered(2, 1);
                if d == org.xy
                    idx = (row-1) * st.xy_n_col + col;
                else
                    idx = (col-1) * st.xy_n_row + row;
                end
                rect = brick.subtract(rect, st.xy_offsets(:, idx));
            end
            
            % x 
            x = sort(rect(1, :));
            x_layout = org.x;
            [x_zoom_dim, x_zoom] = deal([]); % default values
            for ix = length(x_layout):-1:1
                d = x_layout(ix);
                x = x - st.x_offset(ix);
                idx = x / st.x_step(ix);
                
                % if we reached theinternal dimension, things are pretty
                % simple
                if ix == 1
                    x_zoom_dim = d;
                    x_zoom = idx;
                    break
                end
                
                % we accept zooming in an external dimension only if the
                % rectangle covers a full column
                column_width = st.x_span(ix - 1);
                x_in = x - round(idx) * st.x_step(ix);
                idx_inside = round(idx) ...
                    + [(x_in(1) > -column_width/2), -(x_in(2) < column_width/2)];
                if diff(idx_inside) >= 0
                    x_zoom_dim = d;
                    x_zoom = idx_inside + [-.4999 .4999];
                    break
                end
                
                % no valid zoom; we will consider zoom in more internal
                % dimensions only if the y-selection covers, even
                % partially, no more than one row
                idx_covered = round(idx) ...
                    + [(x_in(1) > -column_width/2), -(x_in(2) < column_width/2)];
                if diff(idx_covered)
                    break
                end
                
                % continue: get y-selection coordinates in the next, more
                % internal, dimension
                x = x - round(idx_covered(1)) * st.x_step(ix);
            end
            
            % y 
            y = sort(rect(2, :), 'descend');
            y_layout = org.y;
            [y_zoom_dim, y_zoom] = deal([]); % default values
            broke = false;
            for iy = length(y_layout):-1:1
                d = y_layout(iy);
                y = y - st.y_offset(iy);
                y_step = st.y_step(iy);
                idx = y / y_step;
                
                % if we reached the internal dimension, things are pretty
                % simple
                if iy == 1 && strcmp(G.D.display_mode, 'image')
                    y_zoom_dim = d;
                    y_zoom = idx;
                    broke = true;
                    break
                end
                
                % we accept zooming in an external dimension only if the
                % rectangle covers a full row
                if iy == 1
                    row_height = st.y_available;
                else
                    row_height = st.y_span(iy - 1);
                end
                y_in = y - round(idx) * y_step; % 
                idx_inside = round(idx) ...
                    + [(y_in(1) < row_height/2), -(y_in(2) > -row_height/2)];
                if diff(idx_inside) >= 0
                    y_zoom_dim = d;
                    if isempty(d)
                        y_zoom = [];
                    else
                        y_zoom = idx_inside + [-.4999 .4999];
                    end
                    broke = true;
                    break
                end
                
                % no valid zoom; we will consider zoom in more internal
                % dimensions only if the y-selection covers, even
                % partially, no more than one row
                idx_covered = round(idx) ...
                    + [(y_in(1) < -row_height/2), -(y_in(2) > row_height/2)];
                if diff(idx_covered)
                    broke = true;
                    break
                end
                
                % continue: get y-selection coordinates in the next, more
                % internal, dimension
                y = y - round(idx_covered(1)) * y_step;
            end
            
            % We might be zooming in time courses value as well!
            if ~broke && strcmp(G.D.display_mode, 'time courses')
                % this zoom value will be applied to the current clipping
                clip_zoom = sort(y) / st.y_available + .5; % between 0 and 1 instead of centered on 0
                
                % get also the coordinates in the zslice to know to which
                % display element this clip zoom applies
                zijk = round(G.graph_to_zslice(rect(:,1)));
                clip_zijk = ones(G.D.nd, 1);
                clip_zijk(G.D.clip_dim) = zijk(G.D.clip_dim);
                
                % if we zoom in values, do not zoom in x external
                % dimensions, as this is too disturbing...
                if ~isempty(x_zoom_dim) && x_zoom_dim ~= org.x(1)
                    [x_zoom_dim, x_zoom] = deal([])
                end
            end
                
            % Zooming should occur either in external dimensions, either in
            % internal dimensions or values; not both. If we have both
            % kinds, prefer the internal dimensions.
            zoom_dim = [x_zoom_dim, y_zoom_dim];
            zslice_zoom = [x_zoom; y_zoom];
        end
    end
end

%---
function [x_pair, y_pair] = checkpairs(x_head, y_head, do_signal)
% this function is very ad-hoc and could be improved

nx = length(x_head);
ny = length(y_head);
x_pair = zeros(1, nx);
y_pair = zeros(1, ny);

% restrict to dimensions which are measures
xok = false(1, nx);
yok = false(1, ny);
for i=1:nx, xok(i) = x_head(i).is_measure; end
for i=1:ny, yok(i) = y_head(i).is_measure; end
if ~any(xok) || ~any(yok), return, end

% look for dimensions being in the same space!
xunits = {x_head.unit};
yunits = {y_head.unit};
if xok(1) && yok(1) && ~isempty(xunits{1}) && isequal(xunits(1), yunits(1))
    [x_pair(1), y_pair(1)] = deal(1, 1);
    [xok(1), yok(1)] = deal(false);
end
if nx >= 2 && ny >= 2 && xok(nx) && yok(ny) && ~isempty(xunits{nx}) && isequal(xunits(nx), yunits(ny))
    [x_pair(nx), y_pair(ny)] = deal(ny, nx);
    [xok(nx), yok(ny)] = deal(false);
end
if do_signal && nx >= 2 && xok(2) && yok(1) && ~isempty(xunits{2}) && isequal(xunits(2), yunits(1))
    % (x1 cannot be paired with y2 for images)
    [x_pair(2), y_pair(1)] = deal(1, 2);
    [xok(1), yok(1)] = deal(false);
end
if nx >= 2 && xok(2) && ny >= 2 && yok(2) && ~isempty(xunits{2}) && isequal(xunits(2), yunits(2))
    [x_pair(2), y_pair(2)] = deal(2, 2);
    [xok(2), yok(2)] = deal(false); %#ok<NASGU>
end

end
