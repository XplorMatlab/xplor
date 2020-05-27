classdef ViewDisplay < xplr.GraphNode
    % viewdisplay
    
    % Content
    properties (SetAccess='private')
        % data
        V           % parent 'view' object
        zoom_slicer
        previous_headers = xplr.Header.empty(1, 0);
        % graphics
        hp          % display panel
        ha          % main axes
        menu        % general display menu
        % display
        no_display = false   % data is too large, cancel display
        grid         % containers for line/images that will be translated/scaled
        h_display    % handles of line/images
        grid_clip    % clipping for each grid element
        signals_baseline    % subtracted baseline value for every signal if any
        h_legend     % handle of legend axes
        labels      % xplr.displaylabel object
        graph       % xplr.displaygraph object
        navigation  % xplr.DisplayNavigation object
        clipping    % xplr.ClipTool object
        color_map    % xplr.color_map object
    end
    
    % Some "working memory"
    properties (Access='private')
        slice_change_event
        listeners = struct;
    end
    
    % Display properties
    properties (SetObservable=true, AbortSet=true)
        display_mode = 'image';  % 'time courses' or 'image'
        show_color_legend = false; % false by default because of bug in Matlab's legend function, which inactivates several listeners
        line_alpha = 1; % lines have a small degree of transparency!
        do_title = true;
    end    
    properties (SetAccess='private')
        layout_id                              % layout, i.e. which data dimension appear on which location; set with function set_layout_id
        layout_id_all                           % layout, including singleton dimensions; set with function set_layout_id
        active_dim_id = struct('x', [], 'y', [])   % dimensions on which zooming mouse actions and sliders apply; change with function make_dim_active(D,d)
        color_dim_id = [];                      % set with set_color_dim
    end
    
    % Shortcuts (dependent)
    properties (Dependent, SetAccess='private')
        nd
        slice
        zslice
        zoom_filters
        active_dim
        color_dim
        layout
        internal_dim
        internal_dim_id
        external_dim
        external_dim_id
        external_size
        overlap_dim
        clip_dim
        current_clip
    end
    
    % Constructor, destructor
    methods
        function D = ViewDisplay(V)
            % parent 'view' object and panel
            D.V = V;
            D.hp = V.panels.display;
            set(D.hp, 'deletefcn', @(u,e)delete(V))
            
            % zoom slicer zooms into "slice" to yield "zslice"
            D.zoom_slicer = xplr.ZoomSlicer(V.slicer.slice, D);
            
            % axes
            D.ha = axes('parent', D.hp);
            axis(D.ha, [-.5, .5, -.5, .5]) % center 0, available space 1
            set(D.ha, 'box', 'on', 'clim', [0, 1])
            try set(D.ha, 'XTickLabelRotation', 45), end % recent Matlab versions only
            try set(D.ha, 'TickLabelInterpreter', 'none'), end % recent Matlab versions only
            D.listeners.ax_siz = fn_pixelsizelistener(D.ha, @(u,e)axis_resize(D));
            D.add_listener(D.ha, D.listeners.ax_siz);
            c = disable_listener(D.listeners.ax_siz); % prevent display update following automatic change of axis position during all the following initializations
            
            % 'time courses'/'image' switch
            p = fn_propcontrol(D, 'display_mode', {'popupmenu', 'time courses', 'image'}, 'parent', D.hp);
            fn_controlpositions(p.hu, D.hp, [1, 1], [-90, -25, 90, 25])
            
            % positionning (needed by both labels and data display)
            D.graph = D.add_component(xplr.DisplayGraph(D));
            
            % automatic label positionning
            D.labels = D.add_component(xplr.DisplayLabels(D));
            
            % general display menu
            D.menu = uimenu(D.V.hf, 'label', 'Display', 'callback', @(u,e)display_menu(D));
            
            % clipping tool
            D.clipping = D.add_component(xplr.ClipTool(D)); % creates a menu
            
            % color_map tool
            D.color_map = D.add_component(xplr.ColorMapTool(D)); % creates a menu
            D.add_listener(D.color_map, 'ChangedColorMap', @(u,e)D.update_display('clip'));
            
            % navigation (sliders, mouse actions)
            D.navigation = D.add_component(xplr.DisplayNavigation(D)); % creates a menu
            
            % set organization, connect sliders, display data and labels
            D.slice_change_event = struct('flag', 'global');
            zslice_change(D)
            
            % listeners
            D.add_listener(D.slice, 'ChangedData', @(u,e)set(D, 'slice_change_event', e)); % mark that slice has changed, but treat it only later
            D.add_listener(D.zoom_slicer, 'ChangedZoom', @(u,e)zoom_change(D,e));
            D.add_listener(D.zslice, 'ChangedData', @(u,e)zslice_change(D,e));
            
            % problem: c won't be deleted automatically (and ax_siz listener
            % might not be re-enabled) because the workspace continue to
            % exist, because of all the anonymous functions that were
            % defined
            delete(c)
        end
        function delete(D)
            delete@xplr.GraphNode(D)
            delete(D.zoom_slicer)
            delete(D.navigation)
            delete(D.labels)
            delete(D.graph)
        end
    end
    
    % Dependent properties
    methods
        function n = get.nd(D)
            n = D.V.slicer.slice.nd;
        end
        function x = get.slice(D)
            x = D.V.slicer.slice;
        end
        function x = get.zslice(D)
            x = D.zoom_slicer.slice;
        end
        function zoom_filters = get.zoom_filters(D)
            zoom_filters = [D.zoom_slicer.filters.obj];
        end
        function active_dim = get.active_dim(D)
            active_dim = struct( ...
                'x', D.slice.dimension_number(D.active_dim_id.x), ...
                'y', D.slice.dimension_number(D.active_dim_id.y));
        end
        function color_dim = get.color_dim(D)
            color_dim = D.slice.dimension_number(D.color_dim_id);
            if isequal(color_dim, 0), color_dim = []; end
        end
        function layout = get.layout(D)
            layout = D.layout_id.dimension_number();
        end
        function dim = get.internal_dim(D)
            org = D.layout;
            if isempty(org.x)
                dim = [];
            else
                dim = org.x(1);
            end
            if strcmp(D.display_mode, 'image')
                if ~isempty(org.y)
                    dim = [dim org.y(1), org.merged_data];
                elseif ~isempty(org.merged_data)
                    dim = [dim, org.merged_data];
                end
            end
            dim = sort(dim);
        end
        function dim_id = get.internal_dim_id(D)
            dim = D.internal_dim;
            dim_id = [D.slice.header(dim).dim_id];
        end
        function dim = get.external_dim(D)
            org = D.layout;
            dim = [org.x(2:end), org.y(1+strcmp(D.display_mode,'image'):end) , ...
                org.xy, org.yx];
            dim = sort(dim);
        end
        function dim_id = get.external_dim_id(D)
            dim = D.external_dim;
            dim_id = [D.slice.header(dim).dim_id];
        end
        function sz = get.external_size(D)
            sz = [D.slice.header(D.external_dim).n];
        end
        function dim = get.overlap_dim(D)
            % 'overlap' dimensions are those where several time course
            % signals can appear superimposed within a single grid cell
            if strcmp(D.display_mode,'time courses')
                dim = D.layout.merged_data;
            else
                dim = zeros(1, 0);
            end
        end
        function dim = get.clip_dim(D)
            % mergeddata and external dimensions together: these are the
            % dimensions where clipping range can be defined independently
            org = D.layout;
            dim = [org.x(2:end), org.y(1+strcmp(D.display_mode,'image'):end), ...
                org.xy, org.yx, org.merged_data];
            dim = sort(dim);
        end
    end
    
    % Some usefull simple methods
    methods
        function s = get_size(D, unit, dim)
            % function s = get_size(D, unit [, dim])
            s = fn_objectsize(D.ha, unit);
            if nargin >= 3
                if isnumeric(dim)
                    s = s(dim);
                elseif strcmp(dim, 'x')
                    s = s(1);
                elseif strcmp(dim, 'y')
                    s = s(2);
                else
                    error argument
                end
            end
        end
    end
    
    % Display menu
    methods (Access='private')
        function display_menu(D)
            m = D.menu;
            delete(get(m, 'children'))
            
            % title
            fn_propcontrol(D, 'do_title', 'menu', ...
                {'parent', m, 'label', 'Display title'});

            % time courses display
            if strcmp(D.display_mode, 'time courses')
                fn_propcontrol(D, 'line_alpha', ...
                    {'menuval', {1, .7, .4, .1}, {'none', 'mild', 'medium', 'strong', 'manual'}}, ...
                    {'parent', m, 'label', 'Lines transparency', 'separator', 'on'});
                do_sep = true;
            else
                do_sep = false;
            end
            
            % cross
            fn_propcontrol(D.navigation, 'show_cross', 'menu', ...
                {'parent', m, 'label', 'Show cross', 'separator', onoff(do_sep)});
            if D.navigation.show_cross
                fn_propcontrol(D.navigation,'crosscolor', ...
                    {'menu', {'k', 'b', 'r', [1, 1, 1]*.6, 'w'}, {'black', 'blue', 'red', 'gray', 'white', 'other'}}, ...
                    {'parent',m,'label','Cross color'});
                fn_propcontrol(D.navigation, 'crossalpha', ...
                    {'menu', {1, .4, .05}, {'none', 'medium', 'barely visible', 'manual'}}, ...
                    {'parent', m, 'label', 'Cross transparency'});
            end
            
            % separation marks
            org = D.layout_id;
            if length(org.x) > 1 || length(org.y) > strcmp(D.display_mode, 'image') || ~isempty([org.xy, org.yx])
                fn_propcontrol(D.graph, 'show_separation', 'menu', ...
                    {'parent', m, 'label', 'Show separations between lines/images', 'separator', 'on'});
                if D.graph.show_separation
                    fn_propcontrol(D.graph, 'separationcolor', ...
                        {'menu', {'k', [.8, .8, 1], [1, .8, .8], [.8, .8, .8]}, {'black', 'light blue', 'light red', 'light gray', 'other'}}, ...
                        {'parent', m, 'label', 'Separations color'});
                end
            end            
            
            % when moving dimension label, immediate display update?
            fn_propcontrol(D.labels, 'do_immediate_display', 'menu', ...
                {'parent', m, 'label', 'Immediate display update when moving dimensions', 'separator', 'on'});

            % reset display
            uimenu(m, 'label', 'Reset display', 'separator', 'on', ...
                'callback', @(u,e)D.reset_display())
        end
    end
    
    % Change filters, zoom_filters binning, organization and active dim
    methods (Access='private')
        function out = check_active_dim(D, do_immediate_update, do_auto)
            % function any_chg = check_active_dim(D,do_immediate_update[,do_auto])
            %---
            % check that current active dims are ok, but also set active
            % dims if there aren't
            if nargin<2, do_immediate_update = true; end
            if nargin<3, do_auto = false; end
            
            org_id = D.layout_id;
            sz = D.slice.sz; 
            
            % try to keep same active dims if they still exist in the new
            % data
            dx = D.active_dim_id.x;
            if ~ismember(dx, [D.slice.header.dim_id]), dx = []; end
            dy = D.active_dim_id.y;
            if ~ismember(dy, [D.slice.header.dim_id]), dy = []; end
            
            % check position of active dims
            if isempty([dx dy]) 
                % no valid active dim: set some if do_auto is set to true
                if ~do_auto
                    % do not set active dims
                %                 elseif ~isempty(org_id.xy)
                %                     dy = org_id.xy;
                %                 elseif ~isempty(org_id.yx)
                %                     dx = org_id.yx;
                else
                    % (x)
                    if isempty(dx) && ~isempty(org_id.x)
                        dx = org_id.x(end);
                    end
                    % (y)
                    if isempty(dy) && ~isempty(org_id.y)
                        dy = org_id.y(end);
                    end
                end
            elseif ismember(dx, org_id.x)
                % ok, we can keep dx as the active x dimension
                if ismember(dy, org_id.y)
                    % we can also keep dy as the active y dimension
                    [dx, dy] = deal(dx, dy);
                else
                    [dx, dy] = deal(dx, []);
                end
            elseif ismember(dx, org_id.y)
                % dimension dx moved to y location
                if ismember(dy, org_id.x)
                    % and dimension dy moved to x location: switch the 2
                    % active dims!
                    [dx, dy] = deal(dy, dx);
                elseif ismember(dy, org_id.y)
                    % keep dy as the active y dimension, remove dx from
                    % active x
                    [dx, dy] = deal([], dy);
                else
                    % replace active y dimension with dx
                    [dx, dy] = deal([], dx);
                end
            else
                if ismember(dy, org_id.x)
                    [dx, dy] = deal(dy, []);
                elseif ismember(dy, org_id.y)
                    [dx, dy] = deal([], dy);
                elseif ismember(dx, org_id.xy)
                    [dx, dy] = deal([], dx);
                elseif ismember(dx, org_id.yx)
                    [dx, dy] = deal(dx, []);
                elseif ismember(dy, org_id.xy)
                    [dx, dy] = deal([], dy);
                elseif ismember(dy, org_id.yx)
                    [dx, dy] = deal(dy, []);
                else
                    % should not happen
                    [dx, dy] = deal([], []);
                end
            end
                               
            % update property
            new_value = struct('x', dx, 'y', dy);
            any_chg = ~isequal(new_value, D.active_dim_id);
            if nargout > 0, out = any_chg; end
            if ~any_chg, return, end
            D.active_dim_id = new_value;
            
            % update display
            if do_immediate_update
                D.navigation.connect_zoom_filter()
                D.labels.update_labels('active')
                D.graph.s()
                D.graph.set_value_ticks()
            end
        end
    end    
    methods
        function dimension_context_menu(D, m, dim)
            % function dimension_context_menu(D,m,dim)
            %---
            % This function populates the context menu that appears when
            % right-clicking on a label
            
            [dim, dim_id] = D.slice.dimension_number_and_id(dim);
            head = D.slice.header(dim);
            delete(get(m, 'children'))
            
            % Line properties
            % (color)
            do_color = strcmp(D.display_mode, 'time courses');
            if do_color
                uimenu(m, 'label', ['Color according to ', head.label], 'checked', onoff(isequal(D.color_dim_id,dim_id)), ...
                    'callback', @(u,e)D.set_color_dim(fn_switch(isequal(D.color_dim_id,dim_id),[],dim_id)))
                uimenu(m, 'label', 'Display color legend', ...
                    'enable', onoff(isequal(D.color_dim_id,dim_id)), 'checked', onoff(D.show_color_legend), ...
                    'callback', @(u,e)set(D,'show_color_legend',~D.show_color_legend))
            end

            % Binning
            m1 = uimenu(m, 'label', 'Binning', 'Separator', onoff(do_color));
            bin_values = {1, 2, 3, 4, 5, 10, 20, 'set'};
            bin_displays = {'none', '2', '3', '4', '5', '10', '20', 'other...'};
            cur_bin = D.zoom_filters(dim).bin;
            for i=1:length(bin_values)
                bin = bin_values{i};
                uimenu(m1, 'label', bin_displays{i}, 'checked', onoff(isequal(cur_bin,bin)), ...
                    'callback', @(u,e)set_bin(D,dim,bin));
            end

            % select ZoomFilter key (check the created menu item
            % that corresponds to the current key)
            m2 = uimenu(m, 'label', 'zoom filter', 'Separator', onoff(do_color));
            available_keys = xplr.Bank.available_filter_keys('zoomfilter');
            new_key = max(available_keys) + 1;
            key_values = [0, available_keys, new_key];
            fn_num2str(available_keys, 'shared zoom %i', 'cell');
            key_displays = [ ...
                'private zoom', ...
                fn_num2str(available_keys, 'shared zoom %i', 'cell'), ...
                num2str(new_key,'shared zoom %i (new key)')
                ];
            cur_key = D.zoom_filters(dim).link_key;
            for i=1:length(key_values)
                key_value = key_values(i);
                uimenu(m2, 'label', key_displays{i}, 'checked', onoff(isequal(cur_key,key_value)), ...
                    'callback', @(u,e)D.zoom_slicer.change_key(dim,key_value));
            end

            % select crossSelector key 
            cur_filt = D.navigation.point_filters{dim};
            if ~isempty(cur_filt)
                m2 = uimenu(m, 'label', 'cross selector key', 'Separator', onoff(do_color));

                available_keys = xplr.Bank.available_filter_keys('point');
                new_key = max(available_keys) + 1;
                key_values = [0, available_keys, new_key];
                fn_num2str(available_keys, 'cross selector key %i', 'cell');
                key_displays = [ ...
                    'private cross selector', ...
                    fn_num2str(available_keys, 'cross selector key %i', 'cell'), ...
                    num2str(new_key,'cross selector key %i (new key)')
                    ];
                    cur_key = cur_filt.link_key;
                uimenu(m2, 'label', 'show point selector','callback',@(u,e)xplr.Bank.show_list(cur_filt));
                for i=1:length(key_values)
                    key_value = key_values(i);
                    uimenu(m2, 'label', key_displays{i}, 'checked', onoff(isequal(cur_key,key_value)), ...
                        'callback', @(u,e)connect_point_filter(D.navigation, dim, key_value));
                end
            end
        end
        function set_bin(D,d,bin)
            if strcmp(bin, 'set')
                bin = fn_input('Binning', D.zoom_filters(d).bin, 'stepper 1 1 Inf 1');
                if isempty(bin), return, end
            end
            D.zoom_filters(d).set_bin(bin)
        end
        function set_layout_id(D, new_layout_id, do_immediate_display)
            % function set_layout_id(D,new_layout_id[,do_immediate_display])
            %---
            % if do_immediate_display is set to false, only labels are
            % updated; if it is set to true, update happens regardless of
            % whether newlayout is actually new or not (this allows finishing
            % a previous incomplete update with do_immediate_display set to
            % false)
            if isequal(new_layout_id, D.layout_id_all), return, end
            c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
            if nargin < 3
                do_immediate_display = true;
            end
            D.layout_id_all = new_layout_id;
            D.layout_id = new_layout_id.current_layout(); % keep only dimensions actually displayed
            % is zslice too large for being displayed
            D.check_zslice_size()
            % first update graph (new positionning will be needed for both
            % labels and data display)
            D.graph.compute_steps()
            % update labels
            if do_immediate_display, D.check_active_dim(false), end
            D.labels.update_labels()
            drawnow
            % check whether color dim and active dim remain valid
            D.check_color_dim(false)
            D.check_active_dim(false, true)
            % update ticks and display
            if ~do_immediate_display, return, end
            D.graph.set_ticks()
            update_display(D) % will call set_value_ticks if necessary
            % update slider connections
            connect_zoom_filter(D.navigation)
            % reposition cross
            D.navigation.reposition_cross()
            % update selection display
            D.navigation.display_selection()
        end
        function set_dim_location(D, dim_id, location, do_immediate_display)
            % function set_dim_location(D,dim_id,location,do_immediate_display)
            %---
            % set new location of specific dimensions; locations of other
            % dimensions will automatically be adjusted
            % more details in the help of xplr.DisplayLayout.set_dim_location
            %
            % See also xplr.DisplayLayout.set_dim_location
            
            % new layout
            new_layout_id = D.layout_id_all.set_dim_location(dim_id, location);
            
            % update display
            if nargin < 4, do_immediate_display = true; end
            D.set_layout_id(new_layout_id, do_immediate_display)
        end
        function set.active_dim_id(D,value)
            D.active_dim_id = value;
        end
        function make_dim_active(D, dim_id, flag)
            dim_id = D.slice.dimension_id(dim_id);
            c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
            do_toggle = nargin >= 3 && strcmp(flag, 'toggle');
            % update active dim and connect slider
            if ismember(dim_id, [D.layout_id.x, D.layout_id.yx])
                if do_toggle && any(dim_id == D.active_dim_id.x)
                    D.active_dim_id.x = [];
                else
                    D.active_dim_id.x = dim_id;
                    if ismember(dim_id, D.layout_id.yx) || any(ismember(D.active_dim_id.y, [D.layout_id.xy, D.layout_id.yx]))
                        D.active_dim_id.y = [];
                        D.navigation.connect_zoom_filter('y')
                    end
                end
                D.navigation.connect_zoom_filter('x')
            elseif ismember(dim_id, [D.layout_id.y, D.layout_id.xy])
                if do_toggle && any(dim_id == D.active_dim_id.y)
                    D.active_dim_id.y = [];
                else
                    D.active_dim_id.y = dim_id;
                    if ismember(dim_id, D.layout_id.xy) || any(ismember(D.active_dim_id.x, [D.layout_id.xy, D.layout_id.yx]))
                        D.active_dim_id.x = [];
                        D.navigation.connect_zoom_filter('x')
                    end
                end
                D.navigation.connect_zoom_filter('y')
            end
            % update ticks and labels
            D.graph.set_ticks()
            D.graph.set_value_ticks()
            D.labels.update_labels('active');
        end
    end
    
    % Color
    methods
        function set_color_dim(D, dim, do_immediate_display)
            if ~isnumeric(dim) || numel(dim) > 1, error 'color_dim must be empty or scalar', end
            dim_id = D.slice.dimension_id(dim);
            if isequal(dim_id, D.color_dim_id), return, end
            if ~isempty(dim_id) && ~isempty(D.layout_id.x) && dim_id == D.layout_id.x(1), disp 'first x-dimension cannot be used for colors', return, end
            if nargin < 3, do_immediate_display = true; end
            % set property
            D.color_dim_id = dim_id;
            % update color legend?
            if D.show_color_legend
                display_color_legend(D)
            end
            % update display
            if do_immediate_display && strcmp(D.display_mode, 'time courses')
                D.update_display('color')
            end
        end
        function check_color_dim(D, do_immediate_display)
            c_dim_id = D.color_dim_id;
            if ~isempty(c_dim_id) && ~isempty(D.layout_id.x) && c_dim_id == D.layout_id.x(1)
                % cannot color according to the first x-dimension
                if nargin < 2, do_immediate_display = false; end
                D.set_color_dim([], do_immediate_display)
            end
        end
        function set.line_alpha(D, line_alpha)
            % check
            if ~isscalar(line_alpha) || line_alpha <= 0 || line_alpha > 1, error 'incorrect alpha value for lines', end
            % set value
            D.line_alpha = line_alpha;
            % update display
            if strcmp(D.display_mode, 'time courses')
                D.update_display('color')
            end
        end
        function set.show_color_legend(D, val)
            D.show_color_legend = val;
            display_color_legend(D)
        end
        function display_color_legend(D)
            % delete current legend
            delete(D.h_legend)
            D.h_legend = [];
            % do really display a legend
            d = D.color_dim;
            if isempty(d) || ~D.show_color_legend || strcmp(D.display_mode, 'image'), return, end
            % get line handles
            s = substruct('()', num2cell(ones(1, D.nd)));
            s.subs{d} = ':';
            hl = subsref(D.h_display, s);
            % display legend
            names = D.slice.header(d).get_item_names;
            D.h_legend = fn_colorlegend(row(hl), names, 'SouthWest', 'frame');
        end
    end
    
    % Clipping
    methods
        function clip = get.current_clip(D)
            % get clipping range of the current grid cell
            ijk = D.navigation.get_point_index_position('clip', 'round');
            ijk = row(round(D.graph.slice_to_zslice(ijk)));
            clip = subsref_dim(D.grid_clip, 1+D.clip_dim, ijk(D.clip_dim));
            % center if 'adjust by the mean' mode
            if strcmp(D.display_mode,'time courses') && ~isempty(D.clipping .align_signals)
                clip = clip - feval(D.clipping.align_signals, clip, 1);
            end
        end
        function set_clip(D, clip, do_update_display)
            % input
            if ~isnumeric(clip) || length(clip) ~= 2 || diff(clip) <= 0 || any(isnan(clip)|isinf(clip))
                xplr.debug_info('stop','clip value is not valid')
                return
            end
            if nargin < 3, all_cells = false; end
            % set property
            if all_cells
                D.grid_clip(1,:) = clip(1);
                D.grid_clip(2,:) = clip(2);
            else
                ijk = D.navigation.get_point_index_position('clip',  'round');
                ijk = row(round(D.graph.slice_to_zslice(ijk)));
                dim_indp = D.clipping.independent_dim;
                D.grid_clip = subsasgn_dim(D.grid_clip, [1, 1+dim_indp],[1, ijk(dim_indp)], clip(1));
                D.grid_clip = subsasgn_dim(D.grid_clip, [1, 1+dim_indp],[2, ijk(dim_indp)], clip(2));
            end
            % update display; note that this might also modify D.gridclip
            % if D.cliping.align_signals is not empty
            update_display(D,'clip')
        end
        function auto_clip(D, all_cells)
            if nargin<2, all_cells = false; end
            if all_cells
                D.grid_clip(:) = NaN;
            else
                ijk = D.navigation.get_point_index_position('clip', 'round');
                ijk = row(round(D.graph.slice_to_zslice(ijk)));
                dim_indp = D.clipping.independent_dim;
                D.grid_clip = subsasgn_dim(D.grid_clip, 1+dim_indp, ijk(dim_indp), NaN);
            end
            update_display(D, 'clip')
        end
        function clip = get_clip_range(D,data)
            clip = fn_clip(data(:), D.clipping.auto_clip_mode, 'getrange');
            if isinf(clip(1)), clip(1) = -1e6; end
            if isinf(clip(2)), clip(2) = 1e6; end
        end
    end

    % Title
    methods
        function set.do_title(D, value)
            D.do_title = value;
            D.display_title()
        end
        function display_title(D)
            if ~D.do_title
                title(D.ha,'')
                return
            end
            % names for singleton dimensions
            sz = D.slice.sz;
            singleton_dims = find(sz==1);
            n_singleton = length(singleton_dims);
            names = cell(1, n_singleton);
            for i = 1:n_singleton
                d = singleton_dims(i);
                head = D.slice.header(d);
                namei = head.get_item_name(1);
                % add dimension label if name is an integer
                if mod(str2double(namei), 1)==0
                    namei = [head.label ' ' namei];
                end
                names{i} = namei;
            end
            % add data name if any
            if ~isempty(D.V.data.name)
                names = [D.V.data.name, names];
            end
            title(D.ha, {fn_strcat(names,' - '), ''}) % the second '' serves to place the title slightly higher
        end
    end
    
    % Main display: grid of signals or images
    methods (Static)
        function ok = test_displayable(sz, display_mode, layout)
            sz_grid = sz;
            if ~isempty(layout.x), sz_grid(layout.x(1)) = 1; end
            if strcmp(display_mode, 'image') && ~isempty(layout.y), sz_grid(layout.y(1)) = 1; end
            if strcmp(display_mode, 'time courses')
                ok = (prod(sz_grid) <= xplr.Parameters.get('display.NLineMax')) && ...
                    (prod(sz) <= xplr.Parameters.get('display.NLinePointMax'));
            else
                ok = (prod(sz_grid) <= xplr.Parameters.get('display.NImageMax')) && ...
                    (prod(sz) <= xplr.Parameters.get('display.NImagePixelMax'));
            end
        end
    end
    methods (Access = 'private')
        function update_display(D, flag, dim, ind)
            % function updateDisplay(D[,flag,dim,ind])
            %---
            % available flags are: 'clip' 'global' 'chgdata'
            % 'chgdata&blocksize' 'new' 'remove' 'chg' 'perm' 'pos' 'color'
            if nargin < 3, dim = []; end
            
            % Is data too large for being displayed?
            if D.no_display
                % too many grid elements: cancel display!
                delete_valid(D.grid) % this will also delete children D.h_display
                D.grid_clip = [];
                delete(findall(D.ha, 'type', 'text', 'tag', 'xytick'))
                set(D.ha, 'xtick', [], 'ytick', [])
                if isempty(findall(D.ha, 'type', 'text','tag','no_display'))
                    text(0, 0, {'DATA IS TOO LARGE AND CANNOT BE DISPLAYED', 'BIN IT, OR USE FILTERS TO SLICE IT'}, ...
                        'parent', D.ha, 'horizontalalignment', 'center', 'tag', 'no_display')
                end
                return
            end
            hno_display = findall(D.ha, 'type', 'text', 'tag', 'no_display');
            if ~isempty(hno_display)
                % display was canceled last time: we need a global update
                delete(hno_display)
                flag = 'global';
            end
            
            % Show watch
            c = fn_watch(D.V.hf); %#ok<NASGU>
            
            % To really run fast, avoid accessing object properties
            % repeatedly: access them once for all here
            do_time_courses = strcmp(D.display_mode, 'time courses');
            sz = D.zslice.sz;
            org = D.layout; % convert from dim ID to dim numbers
            internal_dim_ = D.internal_dim;
            external_dim_ = D.external_dim;
            overlap_dim_ = D.overlap_dim;
            clip_dim_ = D.clip_dim;
            elements_dim_ = sort([overlap_dim_ external_dim_]);

            % What to do
            if nargin < 2, flag = 'global'; end
            if ~fn_ismemberstr(flag, {'clip', 'global', 'chg_data', 'chg_data&blocksize', 'new', 'remove', 'chg', 'perm', 'pos', 'color'})
                error 'flag not handled'
            end
            do_reset = strcmp(flag, 'global');
            dim_external = ismember(dim, external_dim_);
            do_new = strcmp(flag, 'new');
            do_remove = strcmp(flag, 'remove');
            do_position = ~fn_ismemberstr(flag, {'chg', 'chg_data', 'clip', 'color'});
            do_data_all = fn_ismemberstr(flag, {'clip', 'global', 'chg_data', 'chg_data&blocksize', 'perm', 'color'}); % color is set when updating ydata, but updating ydata is actually not necessary when only color changes...
            do_data_select = fn_ismemberstr(flag, {'new', 'chg'});
            do_data = do_data_all || do_data_select;
            do_data_select_grid = do_data_select && dim_external;
            do_data_select_overlap = do_data_select && ~dim_external;
            if do_data_select, assert(do_data_select_grid || do_data_select_overlap), end
            do_chg_x = strcmp(flag, 'chg_data&blocksize');
            do_color = do_time_courses && ~fn_ismemberstr(flag, {'chg_data', 'chg_data&blocksize', 'clip'});
            
            % Grid size
            grid_size = ones(1, max(D.nd, 2));
            grid_size(external_dim_) = sz(external_dim_);
            h_display_size = ones(1, max(D.nd, 2));
            h_display_size(elements_dim_) = sz(elements_dim_);
            grid_clip_size = ones(1, max(D.nd, 2));
            grid_clip_size(clip_dim_) = sz(clip_dim_);

            % Check that current grid are valid
            prev_size_check = grid_size;
            if ~isempty(dim)
                prev_size_check(dim) = prev_size_check(dim) + (do_remove-do_new) * length(ind);
            end
            if ~isequal(strict_size(D.grid,length(grid_size)),prev_size_check) ...
                    || ~all(ishandle(D.grid(:)))
                [do_reset, do_position, do_data_all] = deal(true);
                do_color = do_time_courses;
                do_data_select = false;
            end
            
            % Prepare color
            if do_color
                c_dim = D.color_dim;
                color_head = D.zslice.header(c_dim);
                if isempty(D.color_dim)
                    if ~do_reset, set(D.h_display(:), 'color', [0, 0, 0, D.line_alpha]), end
                    do_color = false;
                else
                    k_color = strcmp({color_head.sub_labels.label}, 'ViewColor');
                    if any(k_color)
                        c_map = cell2mat(color_head.values(:,k_color));
                    else
                        c_map = fn_colorset('plot12', 1:color_head.n);
                    end
                    if size(c_map, 2) == 3 && D.line_alpha < 1
                        c_map(:, 4) = D.line_alpha;
                    end
                end
            end
            
            % Prepare display and grid
            if do_reset          % reset display and grid elements
                delete_valid(D.grid) % this will also delete children D.hdisplay
                D.grid = gobjects(grid_size);
                D.grid_clip = NaN([2 grid_clip_size]);
                D.h_display = gobjects(h_display_size);
                [do_position, do_data_all] = deal(true);
            elseif do_new      	% new grid elements
                if ismember(dim,external_dim)
                    D.grid = subsasgn_dim(D.grid, dim, ind, gobjects);
                end
                D.grid_clip = subsasgn_dim(D.grid_clip, 1+dim, ind, NaN);
                D.h_display = subsasgn_dim(D.h_display, dim, ind, gobjects);
            elseif do_remove     % remove grid elements
                if ismember(dim, external_dim)
                    delete_valid(subsref_dim(D.grid, dim, ind)) % this also deletes the children hdisplay objects
                    D.grid = subsasgn_dim(D.grid, dim, ind, []);
                end
                D.grid_clip = subsasgn_dim(D.grid_clip, 1+dim, ind, []);
                delete_valid(subsref_dim(D.h_display, dim, ind)) % this also deletes the children hdisplay objects
                D.h_display = subsasgn_dim(D.h_display, dim, ind, []);
            elseif strcmp(flag, 'chg_data&blocksize') ...
                    && ~do_time_courses && ismember(dim, org.merged_data)
                % D.gridclip might need to be changed in dimension
                % org.mergeddata for images
                % TODO: we have too many cases with this merged_data, need
                % to simplify somehow...
                n_prev = min(size(D.grid_clip, 1+dim), 4);
                n_cur = min(sz(dim), 4);
                if n_cur > n_prev
                    D.grid_clip = subsasgn_dim(D.grid_clip, 1+dim, n_prev+1:n_cur, NaN);
                elseif n_cur < n_prev
                    D.grid_clip = subsasgn_dim(D.grid_clip, 1+dim, n_cur+1:n_prev, []);
                end
            end
            
            % Update clipping: if clipping is independent for each element,
            % postpone clip computation to when these elements will be
            % extracted and displayed. Otherwise compute now.
            if do_data
                displayed_data = D.zslice.data;
                % correct data to align signals on their mean or median?
                align_signals = do_time_courses && D.clipping.align_signals;
                if align_signals
                    D.signals_baseline = displayed_data;
                    if ~isempty(org.x)
                        D.signals_baseline = feval(align_signals, D.signals_baseline, org.x(1));
                        displayed_data = fn_subtract(displayed_data, D.signals_baseline);
                    end
                else
                    D.signals_baseline = [];
                end
                % invalidate clip values that need to be recomputed
                if D.clipping.adjust_to_view && ~strcmp(flag, 'clip')
                    if ismember(dim, D.clipping.independent_dim)
                        D.grid_clip = subsasgn_dim(D.gridclip, 1+dim, ind, NaN);
                    else
                        D.grid_clip(:) = NaN;
                    end
                end
                % compute clip
                dim_indpc = D.clipping.independent_dim; % dimensions with independent clipping
                clip_at_elements_level = isequal(elements_dim_, dim_indpc);
                if strcmp(flag,'perm')
                    % mere permutation
                    D.grid_clip = subsref_dim(D.grid_clip, dim, ind);
                elseif clip_at_elements_level
                    % every element will be clipped independently: postpone
                    % clip computation to when elements will be extracted,
                    % to avoid unneeded repetitive extractions
                else
                    % compute clip independently in sub-parts of the data
                    sz_indpc = ones(1,D.nd);
                    sz_indpc(dim_indpc) = sz(dim_indpc);
                    for idx_indpc = 1:prod(sz_indpc)
                        ijk_indpc = row(fn_indices(sz_indpc(dim_indpc), idx_indpc, 'g2i'));
                        dat_part = subsref_dim(displayed_data, dim_indpc, ijk_indpc);
                        clip_part = subsref_dim(D.grid_clip, 1+dim_indpc, ijk_indpc);
                        idx = find(~isnan(clip_part(1, :)), 1);
                        if isempty(idx)
                            % need to recompute
                            clip = D.get_clip_range(dat_part);
                        else
                            % use current clip value where it is not
                            % defined
                            clip = clip_part(:, idx);
                        end
                        D.grid_clip = subsasgn_dim(D.grid_clip, [1 1+dim_indpc], [1 ijk_indpc], clip(1));
                        D.grid_clip = subsasgn_dim(D.grid_clip, [1 1+dim_indpc], [2 ijk_indpc], clip(2));
                    end
                end
            end
            
            % Create or modify grid containers
            if do_position
                % (list of grid cell indices)
                idx_grid_list = 1:prod(grid_size);
                if do_data_select_grid && ~doposition
                    % not all grid elements need to be visited ('chg' flag)
                    idx_grid_list = reshape(idx_grid_list, grid_size);
                    subs = substruct('()', repmat({':'}, 1, D.zslice.nd));
                    subs.subs{dim} = ind;
                    idx_grid_list = row(subsref(idx_grid_list, subs));
                end
                ijk_grid_list = fn_indices(grid_size, idx_grid_list, 'g2i');
                % (prepare dispatch)
                M = D.graph.get_transform(ijk_grid_list);
                % (loop on grid cells)
                for u = 1:length(idx_grid_list)
                    idx_grid = idx_grid_list(u);
                    ijk_grid = ijk_grid_list(:,u);
                    do_data_cell = do_data && (~do_data_select_grid || any (ind == ijk_grid(dim)));
                    do_create_cur = do_reset || (do_new && do_data_cell);
                    if do_create_cur
                        D.grid(idx_grid) = hgtransform('parent', D.ha, ...
                            'matrix', M(:, :, idx_grid), 'HitTest', 'off');
                    else
                        set(D.grid(idx_grid), 'matrix', M(:, :, idx_grid))
                    end
                end
            end

            % Create or modify individual signals/images
            if do_data
                % list of all signals/images indices
                idx_h_display_list = reshape(1:prod(h_display_size), h_display_size);
                if do_data_select
                    % not all objects need to be visited
                    subs = substruct('()', repmat({':'}, 1, D.nd));
                    subs.subs{dim} = ind;
                    idx_h_display_list = subsref(idx_h_display_list, subs);
                end
                % (reorganize list as: rows=different grid cells, columns=overlapped elements inside same grid cell)
                idx_h_display_list = fn_reshapepermute(idx_h_display_list, {external_dim_ overlap_dim_ internal_dim_});
                [n_grid, n_overlap] = size(idx_h_display_list); % can be smaller than total numbers of grids/overlaps
                ijk_h_display_list = fn_indices(h_display_size, idx_h_display_list(:),'g2i');
                ijk_h_display_list = reshape(ijk_h_display_list, [D.nd, n_grid, n_overlap]);
                % (corresponding grid cells)
                ijk_grid_list = ijk_h_display_list(:, :, 1);
                ijk_grid_list(overlap_dim_, :) = 1;
                idx_grid_list = fn_indices(grid_size, ijk_grid_list, 'i2g');

                % subs structure for slicing
                subs = substruct('()', repmat({':'}, 1, length(sz)));
                % (and permutation that follows)
                nd_perm = max(3, D.nd);
                internal_perm = zeros(1, nd_perm);
                if do_time_courses
                    % reorder dimensions as x-others
                    if ~isempty(org.x), internal_perm(1) = org.x(1); end
                else
                    % reorder dimensions as y-x-channel-others
                    if ~isempty(org.y), internal_perm(1) = org.y(1); end
                    if ~isempty(org.x), internal_perm(2) = org.x(1); end
                    if ~isempty(org.merged_data), internal_perm(3) = org .merged_data; end
                end
                internal_perm(~internal_perm) = setdiff(1:nd_perm, internal_dim_);

                % go! loop on grid cells and on elements inside these cells
                for u = 1:n_grid
                    idx_grid = idx_grid_list(u);
                    for v = 1:n_overlap
                        idx_h_display = idx_h_display_list(u, v);
                        ijk_h_display = ijk_h_display_list(:, u, v);

                        % get the data and compute clipping if necessary
                        [subs.subs{elements_dim_}] = dealc(ijk_h_display (elements_dim_));
                        xi = subsref(displayed_data, subs);
                        xi = permute(xi, internal_perm);
                        xi = fn_float(xi);
                        clipi = subsref_dim(D.grid_clip, 1+elements_dim_, ...
                            ijk_h_display(elements_dim_)); % 2x1 vector, or 2xn array if color image and color channel clipped independently
                        if clip_at_elements_level && any(isnan(clipi(:)))
                            % independent clip for this element
                            clipi = D.get_clip_range(xi);
                            clipi = repmat(clipi(:), [1, size(xi,3)]);
                            D.grid_clip = subsasgn_dim(D.grid_clip, 1+elements_dim_, ...
                                ijk_h_display(elements_dim_), clipi);
                        end
                        % display it
                        if do_time_courses
                            xi = (xi - clipi(1)) / diff(clipi);
                            nt = size(xi,1);
                            if nt==1
                                line_opt = {'linestyle', 'none', 'marker', ' .'};
                            else
                                line_opt = {'linestyle', '-', 'marker', 'none'};
                            end
                            if ~ishandle(D.h_display(idx_hdisplay))
                                hl = line(1:nt,xi, ...
                                    'parent', D.grid(idx_grid), 'HitTest', 'off', line_opt{:});
                                D.h_display(idx_h_display) = hl;
                            else
                                hl = D.h_display(idx_h_display);
                                if do_chg_x
                                    set(hl, 'xdata', 1:nt, line_opt{:})
                                end
                                set(hl, 'ydata', xi)
                            end
                            if do_color
                                set(hl, 'color', cmap(ijk_h_display(cdim), :))
                            end
                        else
                            % size in color dimension must be 1, 3 or 4;
                            % correct if it is not the case
                            alpha = [];
                            nc = size(xi, 3);
                            if nc > 4
                                xi = xi(:, :, 1:4);
                            end
                            im = D.color_map.color_image(xi, clipi);
                            if nc == 2
                                % add a third blue channel set to zero
                                im(:, :, 3) = 0;
                            elseif nc == 4
                                [im, alpha] = deal(im(:, :, 1:3), im(:, :, 4));
                            end
                            if ~ishandle(D.h_display(idx_h_display))
                                % y coordinates are negative to orient the
                                % image downward (see also comment inside of
                                % displaygaph.gettransform method, where the
                                % y-scale of the hgtransform cannot be set
                                % negative)
                                D.h_display(idx_h_display) = surface( ...
                                    [.5, size(im,2)+.5], [-.5, -.5-size(im,1)], zeros(2), ...
                                    'parent', D.grid(idx_grid), ...
                                    'EdgeColor', 'none', 'FaceColor', 'texturemap', ...
                                    'CDataMapping', 'scaled', 'FaceAlpha', 'texturemap', ...
                                    'CData', im, 'AlphaData', alpha, ...
                                    'HitTest', 'off');
                            elseif do_chg_x
                                set(D.h_display(idx_h_display), ...
                                    'xdata', [.5, size(im,2)+.5], 'ydata', [-.5, -.5-size(im, 1)], ...
                                    'CData', im, 'AlphaData', alpha)
                            else
                                set(D.h_display(idx_h_display), 'CData', im, 'AlphaData', alpha)
                            end
                        end
                    end
                end
            end
            
            % update value y-ticks
            if do_data
                D.graph.set_value_ticks()
            end
            
            % make sur containers are below labels, selections, etc.
            % (fater to have only one call to uistack rather than after
            % creating each element)
            if do_reset || do_new
                uistack(D.grid(:), 'bottom')
            end

        end
        function check_zslice_size(D)
            D.no_display = ~D.test_displayable(D.zslice.sz, D.display_mode, D.layout);
        end
    end
    methods
        function reset_display(D)
            % reset axis
            cla(D.ha)
            % re-display everything
            D.slice_change_event = struct('flag', 'global');
            D.navigation.display_selection('reset')
            zslice_change(D) % this will automatically re-create the cross, but not the selection displays
        end
    end
    
    % Main routines: response to changes in slice, zslice, zoom, etc.
    methods
        function zslice_change(D, e)
            if nargin < 2, flag = 'global'; else flag = e.flag; end
            xplr.debug_info('viewdisplay', 'zslice_change %s', flag)
            c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
            
            % Did slice change as well?
            if ~isempty(D.slice_change_event)
                slice_change(D,D.slice_change_event)
                D.slice_change_event = [];
            end
            
            % Dimension(s) where change occured
            if nargin >= 2
                [chg_dim, chg_dim_id] = D.zslice.dimension_number_and_id(e.dim);
            end
            
            % Is zslice too large for being displayed
            D.check_zslice_size()
            
            % Update graph (will be needed by both labels and data display)
            prevsz = D.graph.zslice_sz;
            D.graph.compute_steps()
            
            % Update labels and ticks (do the labels first because ticks
            % update can change the size of the axes, and therefore trigger
            % labels re-positionning, which can cause error if the number
            % of labels has decreased)
            if fn_ismemberstr(flag, {'all', 'new', 'remove', 'chg&new', 'chg&rm', 'global', 'chg_dim'})
                switch flag
                    case 'global'
                        D.labels.update_labels('global')
                    case 'chg_dim'
                        D.labels.update_labels(flag, chg_dim)
                    otherwise
                        D.labels.update_labels()
                end
            end
            if ~(strcmp(flag,'chg_data') || (strcmp(flag, 'chg') && ~any(chg_dim == [D.active_dim.x, D.active_dim.y])))
                D.graph.set_ticks()
            end

            % Update display
            if fn_ismemberstr(flag, {'global', 'chg_dim'})
                % Reset display
                update_display(D, 'global')
            elseif strcmp(flag, 'chg_data')
                % No change in size, all data need to be redisplayed
                update_display(D, 'chg_data')
            else
                % Smart display update
                if ismember(chg_dim_id, D.internal_dim_id)
                    % changes are within elements (the grid arrangement
                    % remains the same)
                    if fn_ismemberstr(flag, {'perm', 'chg'}) ...
                            || (strcmp(flag, 'all') && D.zslice.header(chg_dim).n == prevsz(chg_dim))
                        flag = 'chg_data'; % no change in size
                    else
                        flag = 'chg_data&blocksize';
                    end
                    D.update_display(flag, chg_dim, e.ind)
                else
                    % the grid arrangement changes
                    switch flag
                        case {'chg', 'new', 'remove', 'perm'}
                            update_display(D, flag, chg_dim, e.ind)
                        case 'all'
                            n_cur = size(D.grid, chg_dim);
                            n = D.zslice.sz(chg_dim);
                            if n == n_cur
                                update_display(D, 'chg_data')
                            elseif n > n_cur
                                update_display(D, 'new', chg_dim, n_cur+1:n)
                                update_display(D, 'chg', chg_dim, 1:n_cur)
                            else
                                update_display(D, 'remove', chg_dim, n+1:n_cur)
                                update_display(D, 'chg', chg_dim, 1:n)
                            end
                        case 'chg&new'
                            update_display(D, 'new', chg_dim, e.ind{2})
                            update_display(D, 'chg', chg_dim, e.ind{1})
                        case 'chg&rm'
                            update_display(D, 'remove', chg_dim, e.ind{2})
                            update_display(D, 'chg', chg_dim, e.ind{1})
                        otherwise
                            error('flag ''%s'' is not handled', 'flag')
                    end
                end
            end
                      
           % reposition cross
           D.navigation.reposition_cross()
           
           % update selections display
           D.navigation.display_selection('referentialchanged')

            % Update legend
            if strcmp(flag, 'global')
                D.color_dim_id = [];
            elseif ~strcmp(flag, 'chg_data') && isequal(chg_dim, D.color_dim)
                display_color_legend(D)
            end
        end
        function slice_change(D, e)
            % function slice_change(D, e)
            %---
            % Function slicechange takes care of changes specific to a
            % slice change (e.g. updating the layout, slider connections).
            % It does not take care of changes common with a zslice change,
            % because it is called from within the zslicechange method,
            % which will take care or them afterwards. Note that the
            % zslicechange method itself was called following the
            % zoomslicer update that was triggered by the slice change,
            % which means that the zoomslicer is already set correctly.

            if nargin < 2, flag = 'global'; else, flag = e.flag; end
            xplr.debug_info('viewdisplay', 'slice_change %s', flag)

            % first time?
            if isempty(D.layout_id_all)
                % some heuristics to choose initial layout
                D.display_mode = fn_switch(sum(D.slice.sz>1) == 1, 'time courses', 'image');
                D.layout_id_all = xplr.DisplayLayout(D);
                D.layout_id = D.layout_id_all;
            else
                % keep locations of dimensions already present in
                % D.layout_id_all, use some heuristic to choose
                % locations of new dimensions
                [D.layout_id_all, D.layout_id] = D.layout_id_all.update_layout();
            end

            % Update active dim and slider connections
            if fn_ismemberstr(flag, {'global'})
                D.check_active_dim(false, true)
                D.navigation.connect_zoom_filter()
            elseif fn_ismemberstr(flag, {'chg_data', 'chg'})
                % slice size did not change
            else
                D.check_active_dim(false)
                D.navigation.connect_zoom_filter()
            end

            % Update color dim
            D.check_color_dim(false)

            % Assign point filters to each updated dimension
            switch flag
                case 'chg_data'
                    % nothing to do: only the data has changed
                case 'global'
                    D.navigation.connect_point_filter()
                case {'all', 'chg_dim', 'new', 'remove', 'chg', 'chg&new', 'chg&rm', 'perm'}
                    dim = D.slice.dimension_number(e.dim);
                    D.navigation.connect_point_filter(dim)
                otherwise
                    error('flag ''%s'' not handled', flag)
            end

            % Check whether current dimension for selections display is
            % still valid, i.e. whether the connected filter still fits the
            % dimension in the new slice (if it is still valid, note that
            % selection display update will occur in D.zslice_change)
            D.navigation.check_selection_filter()

            % Se previous headers to current headers
            D.previous_headers = D.slice.header;

            % Update title
            D.display_title()
        end
        function zoom_change(D, e)
            % update graph positions: if data has changed in size,
            % positioning will be updated upon notification of data change;
            % however if data has not changed in size, positioning needs to
            % be updated here
            if ~e.chg_n_out
                c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
                D.check_zslice_size % is zslice too large for being displayed
                D.graph.compute_steps()
                D.graph.set_ticks()
                D.graph.set_value_ticks()
                D.labels.update_labels()
                update_display(D, 'pos')
                D.navigation.reposition_cross()
                D.navigation.display_selection()
            end
        end
        function axis_resize(D)
            c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
            D.graph.compute_steps()
            D.graph.set_ticks()
            D.graph.set_value_ticks()
            update_display(D, 'pos')
            D.navigation.reposition_cross()
            D.navigation.display_selection()
        end
        function set.display_mode(D, mode)
            c = disable_listener(D.listeners.ax_siz); %#ok<MCSUP,NASGU> % prevent display update following automatic change of axis position
            % set property
            D.display_mode = mode;
            % for 'image' mode, check that layout is valid, and modify it if
            % necessary
            if strcmp(mode, 'image') && length(D.layout_id.merged_data) > 1
                % it is not possible to have more than 1 color dimension ->
                % move other dimensions at 'merged_data' location to 'y'
                new_layout_id = D.layout_idmem;
                dim_id_move = D.layout_id.merged_data(2:end);
                new_layout_id.merged_data = setdiff(new_layout_id.merged_data, dim_id_move, 'stable');
                new_layout_id.y = [new_layout_id.y, dim_id_move];
                D.set_layout_id(new_layout_id) % this automatically updates display among other things
            else
                % update display
                D.check_zslice_size() % is zslice too large for being displayed
                D.graph.compute_steps() %#ok<MCSUP>
                D.graph.set_ticks() %#ok<MCSUP>
                D.labels.update_labels() %#ok<MCSUP>
                update_display(D) % will call set_value_ticks
                D.navigation.display_selection() %#ok<MCSUP>
            end
            % show/hide color legend
            display_color_legend(D)
        end
    end
    
end
