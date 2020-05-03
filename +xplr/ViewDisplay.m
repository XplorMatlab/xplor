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
        h_transform  % containers for line/images that will be translated/scaled
        h_display    % handles of line/images
        grid_clip    % clipping for each grid element
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
    end    
    properties (SetAccess='private')
        layout_id                              % layout, i.e. which data dimension appear on which location; set with function set_layout_id
        layout_id_all                           % layout, including singleton dimensions; set with function set_layout_id
        active_dim_id = struct('x', [], 'y', [])   % dimensions on which zooming mouse actions and sliders apply; change with function make_dim_active(D,d)
        color_dim_id = [];                      % set with set_color_dim
        clip = [0, 1]                          % set with set_clip, auto-clip with auto_clip, other clip settings with sub-object cliptool
    end
    
    % Shortcuts (dependent)
    properties (Dependent, SetAccess='private')
        nd
        slice
        z_slice
        zoom_filters
        active_dim
        color_dim
        layout
        internal_dim
        internal_dim_id
    end
    
    % Constructor, destructor
    methods
        function D = ViewDisplay(V)
            % parent 'view' object and panel
            D.V = V;
            D.hp = V.panels.display;
            set(D.hp, 'deletefcn', @(u,e)delete(V))
            
            % zoom slicer zooms into "slice" to yield "z_slice"
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
            D.clipping = D.add_component(xplr.ClipTool(V.hf)); % creates a menu
            D.add_listener(D.clipping, 'ChangedClip', @(u,e)clip_change(D,e));
            
            % color_map tool
            D.color_map = D.add_component(xplr.colorMapTool(D)); % creates a menu
            D.add_listener(D.color_map, 'ChangedColorMap', @(u,e)D.update_display('clip'));
            
            % navigation (sliders, mouse actions)
            D.navigation = D.add_component(xplr.DisplayNavigation(D)); % creates a menu
            
            % set organization, connect sliders, display data and labels
            D.slice_change_event = struct('flag', 'global');
            z_slice_change(D)
            
            % listeners
            D.add_listener(D.slice, 'ChangedData', @(u,e)set(D, 'slice_change_event', e)); % mark that slice has changed, but treat it only later
            D.add_listener(D.zoom_slicer, 'ChangedZoom', @(u,e)zoom_change(D,e));
            D.add_listener(D.z_slice, 'ChangedData', @(u,e)z_slice_change(D,e));
            
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
        function x = get.z_slice(D)
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
            else
            end
        end
        function dim_id = get.internal_dim_id(D)
            dim = D.internal_dim();
            dim_id = [D.slice.header(dim).dim_id];
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
    
    % General methods
    methods (Access='private')
        function display_menu(D)
            m = D.menu;
            delete(get(m, 'children'))
            
            % time courses display
            if strcmp(D.display_mode, 'time courses')
                fn_propcontrol(D, 'line_alpha', ...
                    {'menuval', {1, .7, .4, .1}, {'none', 'mild', 'medium', 'strong', 'manual'}}, ...
                    {'parent', m, 'label', 'Lines transparency'});
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
            c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
            if nargin < 3
                if isequal(new_layout_id, D.layout_idmem), return, end
                do_immediate_display = true;
            end
            D.layout_id_all = new_layout_id;
            D.layout_id = new_layout_id.current_layout(); % keep only dimensions actually displayed
            % is z_slice too large for being displayed
            D.check_z_slice_size()
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
        function set_clip(D, clip, do_update_display)
            if ~isnumeric(clip) || length(clip) ~= 2 || diff(clip) <= 0 || any(isnan(clip)|isinf(clip))
                xplr.debug_info('stop','clip value is not valid')
                return
            end
            if all(clip == D.clip), return, end
            if nargin < 3, do_update_display = true; end
            % set property
            D.clip = double(clip);
            % update display
            if do_update_display, update_display(D, 'clip'), end
        end
        function auto_clip(D, do_update_display)
            if nargin < 2, do_update_display = true; end
            try
                val = fn_clip(D.z_slice.data(:), D.clipping.auto_clip_mode, 'getrange');
                if isinf(val(1)), val(1) = -1e6; end
                if isinf(val(2)), val(2) = 1e6; end
                if ~any(isnan(val)), set_clip(D,val,do_update_display), end
            catch ME
                disp(ME)
            end
        end
        function clip_change(D, e)
            switch e.flag
                case 'clip'
                    D.set_clip(e.value)
                case 'automode'
                    D.auto_clip()
                case 'adjust'
                    D.update_display('clip')
                case 'span'
                    if strcmp(D.clipping.span, 'curview')
                        D.auto_clip()
                    end
                otherwise
                    error('invalid ChangedClip flag ''%s''', e.flag)
            end
        end
    end
    
    % Update display
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
        function slice_change(D, e)
            % function slice_change(D,e)
            %---
            % function slice_change updates the 'layout' property and slider
            % connections, but does not update display
            % it should be called after the zoom_slicer has already
            % completed its update after slice change (hence the use of the
            % 'slice_change_event' property to delay call to slice_change)
            
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
            elseif fn_ismemberstr(flag, {'chgdata', 'chg'})
                % slice size did not change
            else
                D.check_active_dim(false)
                D.navigation.connect_zoom_filter()
            end
            
            % Update color dim
            D.check_color_dim(false)
            
            % Assign point filters to each updated dimension
            switch flag
                case 'chgdata'
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
            % selection display update will occur in D.z_slice_change)
            D.navigation.check_selection_filter()
            
            % Se previous headers to current headers
            D.previous_headers = D.slice.header;
        end
        function update_display(D, flag, dim, ind)
            % function update_display(D[,flag,dim,ind])
            if nargin < 3, dim = []; end
            
            % Is data too large for being displayed?
            if D.no_display
                % too many grid elements: cancel display!
                deleteValid(D.h_transform) % this will also delete children D.h_display
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
            clip_adjust = D.clipping.adjust;
            
            % What to do
            if nargin < 2, flag = 'global'; end
            if ~fn_ismemberstr(flag, {'clip', 'global', 'chgdata', 'chgdata&blocksize', 'new', 'remove', 'chg', 'perm', 'pos', 'color'})
                error 'flag not handled'
            end
            do_reset = strcmp(flag, 'global');
            do_new = strcmp(flag, 'new');
            do_remove = strcmp(flag, 'remove');
            do_position = ~fn_ismemberstr(flag, {'chg', 'chgdata', 'clip', 'color'});
            do_data_all = fn_ismemberstr(flag, {'clip', 'global', 'chgdata', 'chgdata&blocksize', 'perm', 'color'}); % color is set when updating ydata, but updating ydata is actually not necessary when only color changes...
            do_data_select = fn_ismemberstr(flag, {'new', 'chg'});
            do_data = do_data_all || do_data_select;
            do_chg_x = strcmp(flag, 'chgdata&blocksize');
            do_color = do_time_courses && ~fn_ismemberstr(flag, {'chgdata', 'chgdata&blocksize', 'clip'});
            
            % Reshape slice data adequately: after reshape, dimensions in x
            % will alternate between external and internal dimensions
            sz = D.z_slice.sz;
            org = D.layout; % convert from dim ID to dim numbers
            if ~isempty(org.x)
                x_layout_0 = D.layout.x(1);
            else
                % no dimension displayed in x; mimick an additional
                % dimension displayed in x
                x_layout_0 = length(sz) + 1;
                sz(end+1) = 1;
            end
            % build a 3-element vector of 'internal' dimensions; some of
            % them might be fake dimensions, the idea is to have in any
            % case 3 internal dimensions to avoid handling multiple
            % separate cases
            if do_time_courses
                internal_dim = [x_layout_0 length(sz) + (1:2)];
                sz(end+(1:2)) = 1;
            else
                if ~isempty(org.y)
                    y_layout_0 = org.y(1);
                else
                    y_layout_0 = length(sz) + 1;
                    sz(end+1) = 1;
                end
                if ~isempty(org.merged_data)
                    channel_dim = org.merged_data;
                else
                    channel_dim = length(sz)+1;
                    sz(end+1) = 1;
                end
                internal_dim = [x_layout_0, y_layout_0, channel_dim];
            end
            % prepare the dimensions permutation after extracting each
            % sub-signal/sub-image
            [internal_dim_sort, perm_xi_2_slice] = sort(internal_dim);
            perm_slice_2_xi = [0, 0, 0, 1, 3, 5, 7];
            perm_slice_2_xi(perm_xi_2_slice) = [2, 4, 6];
            if ~do_time_courses
                % for image, XPLOR dimension order is x-y, but Matlab is
                % y-x, so we permute the two before displaying images
                perm_slice_2_xi(1:2) = perm_slice_2_xi([2, 1]);
            end
            % size of reshape
            szr = zeros(1, 7);
            for i = 1:3
                if i == 1
                    szr(2*i-1) = prod(sz(1:internal_dim_sort(1)-1));
                else
                    szr(2*i-1) = prod(sz(internal_dim_sort(i-1)+1:internal_dim_sort(i)-1));
                end
                szr(2*i) = sz(internal_dim_sort(i));
            end
            szr(7) = prod(sz(internal_dim_sort(3)+1:end));
            % reshape
            x = reshape(D.z_slice.data, szr);
            % size of remaining external dimensions
            szo = szr(1:2:end);

            % Check that current h_transform and h_display are valid
            nd = D.z_slice.nd;
            sz1 = sz;
            sz1(internal_dim) = 1;
            sz1prev = sz1;
            if ~isempty(dim), sz1prev(dim) = sz1prev(dim) + (do_remove-do_new)*length(ind); end
            if ~isequal(strictsize(D.h_transform,nd),sz1prev) || ~all(ishandle(D.h_transform(:)))...
                    || ~isequal(strictsize(D.h_transform,nd),sz1prev) || ~all(ishandle(D.h_transform(:)))
                [do_reset, do_position, do_data_all] = deal(true);
                do_data_select = false;
            end
            
            % Prepare color
            if do_color
                c_dim = D.color_dim;
                color_head = D.z_slice.header(c_dim);
                if isempty(D.color_dim)
                    if ~do_reset, set(D.h_display(:), 'color', [0, 0, 0, D.line_alpha]), end
                    do_color = false;
                else
                    k_color = strcmp({color_head.sublabels.label}, 'ViewColor');
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
            sz1 = sz;
            sz1(internal_dim) = 1;
            if do_reset          % reset display and grid elements
                deleteValid(D.h_transform) % this will also delete children D.h_display
                [D.h_transform, D.h_display] = deal(g_objects([sz1, 1]));
                D.grid_clip = zeros([2, sz1, 1]);
                [do_position, do_data_all] = deal(true);
            elseif do_new      	% new grid elements
                subs = substruct('()', repmat({':'}, 1, D.z_slice.nd));
                subs.subs{dim} = ind;
                D.h_transform = subsasgn(D.h_transform, subs, g_objects);
                D.h_display = subsasgn(D.h_display, subs, g_objects);
                subs.subs = [{':'}, subs.subs];
                D.grid_clip = subsasgn(D.grid_clip,subs, 0);
            elseif do_remove     % remove grid elements
                subs = substruct('()', repmat({':'}, 1, D.z_slice.nd));
                subs.subs{dim} = ind;
                deleteValid(subsref(D.h_transform, subs)) % this also deletes the children h_display objects
                D.h_transform = subsasgn(D.h_transform,subs,[]);
                D.h_display = subsasgn(D.h_display,subs,[]);
                subs.subs = [{':'} subs.subs];
                D.grid_clip = subsasgn(D.grid_clip,subs,[]);
            end
            
            % Prepare clipping
            clip_0 = D.clip;
            clip_extent = diff(clip_0);
            
            % Prepare several list of indices beforehand to avoid repeated
            % calls to functions such as ind_2_sub
            idx_list = 1:prod(sz1);
            if do_data_select && ~do_position
                % not all grid elements need to be visited ('chg' flag)
                idx_list = reshape(idx_list, sz1);
                subs = substruct('()', repmat({':'}, 1, D.z_slice.nd));
                subs.subs{dim} = ind;
                idx_list = row(subsref(idx_list, subs));
            end
            ijk_list = fn_indices(sz1, idx_list, 'g2i');
            [i1, i2, i3, i4] = ind_2_sub(szo, idx_list);
            
            % Prepare dispatch
            if do_position
                M = D.graph.get_transform(ijk_list);
            end
            
            % Go! Loop on grid elements
            for u = 1:length(idx_list)
                idx = idx_list(u);
                ijk = ijk_list(:, u);
                do_data_cur = do_data_all || (do_data_select && any(ind==ijk(dim)));
                do_create_cur = do_reset || (do_new && do_data_cur);
                % container
                if do_position
                    if do_create_cur
                        D.h_transform(idx) = hgtransform('parent', D.ha, 'matrix', M(:,:,idx), 'HitTest', 'off');
                    else
                        set(D.h_transform(idx), 'matrix', M(:,:,idx))
                    end
                end
                % line/image
                if do_create_cur || do_data_cur
                    % get the data and adjust clipping if requested
                    xi = x(i1(u), :, i2(u), :, i3(u), :, i4(u));
                    xi = permute(xi, perm_slice_2_xi);
                    switch clip_adjust
                        case 'none'
                            clip_i = clip_0;
                        case 'mean'
                            % adjustment by the mean
                            clip_i = nmean(xi(:)) + [-.5 +.5] * clip_extent;
                        otherwise
                            error('unknown clipping adjustment flag ''%s''', clip_adjust)
                    end
                    xi = fn_float(xi);
                    clip_i = clip_i;
                    % store clipping values
                    D.grid_clip(:, idx) = clip_i;
                    % display it
                    if do_time_courses
                        xi = (xi-clip_i(1))/clip_extent;
                        nt = sz(internal_dim(1));
                        if nt == 1
                            line_opt = {'linestyle', 'none', 'marker', '.'};
                        else
                            line_opt = {'linestyle', '-', 'marker', 'none'};
                        end
                        if do_create_cur
                            hl = line(1:nt, xi, ...
                                'parent', D.h_transform(idx), 'HitTest', 'off', line_opt{:});
                            D.h_display(idx) = hl;
                            if do_color, set(hl, 'color', c_map(ijk(c_dim),:)), end
                        else
                            hl = D.h_display(idx);
                            if do_chg_x
                                set(hl, 'xdata', 1:nt,line_opt{:})
                            end
                            set(hl, 'ydata', xi)
                            if do_color, set(hl, 'color', c_map(ijk(c_dim),:)), end
                        end
                    else
                        % size in color dimension must be 1, 3 or 4;
                        % correct if it is not the case
                        alpha = [];
                        nc = sz(internal_dim(3));
                        if nc > 4
                            xi = xi(:, :, 1:3);
                        end
                        im = D.color_map.color_image(xi, clip_i);
                        if nc == 2
                            % add a third blue channel set to zero
                            im(:, :, 3) = 0;
                        elseif nc == 4
                            [im, alpha] = deal(im(:, :, 1:3), im(:, :, 4));
                        end
                        if do_create_cur
                            % y coordinates are negative to orient the
                            % image downward (see also comment inside of
                            % displaygaph.get_transform method, where the
                            % y-scale of the hgtransform cannot be set
                            % negative)
                            D.h_display(idx) = surface([.5, size(im,2)+.5], [-.5, -.5-size(im,1)], zeros(2), ...
                                'parent', D.h_transform(idx), ...
                                'EdgeColor', 'none', 'FaceColor', 'texturemap', ...
                                'CDataMapping', 'scaled', 'FaceAlpha', 'texturemap', ...
                                'CData', im, 'AlphaData', alpha, ...
                                'HitTest', 'off');
                        elseif do_chg_x
                            set(D.h_display(idx), 'xdata', [.5, size(im,2)+.5], 'ydata', [-.5, -.5-size(im,1)], ...
                                'CData', im, 'AlphaData', alpha)
                        else
                            set(D.h_display(idx), 'CData', im, 'AlphaData', alpha)
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
                uistack(D.h_transform(:), 'bottom')
            end

        end
        function check_z_slice_size(D)
            D.no_display = ~D.test_displayable(D.z_slice.sz, D.display_mode, D.layout);
        end
    end
    methods
        function reset_display(D)
            % reset axis
            cla(D.ha)
            % re-display everything
%             D.slice_change_event = struct('flag','global');
            D.navigation.display_selection('reset')
            z_slice_change(D) % this will automatically re-create the cross, but not the selection displays
        end
    end
    
    % Callbacks (i.e. handle specific changes in slice, properties, etc.)
    methods
        function z_slice_change(D, e)
            if nargin < 2, flag = 'global'; else flag = e.flag; end
            xplr.debug_info('viewdisplay', 'z_slice_change %s', flag)
            c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
            
            % Did slice change as well?
            if ~isempty(D.slice_change_event)
                slice_change(D,D.slice_change_event)
                D.slice_change_event = [];
            end
            
            % Dimension(s) where change occured
            if nargin >= 2
                [chg_dim, chg_dim_id] = D.z_slice.dimension_number_and_id(e.dim);
            end
            
            % Is z_slice too large for being displayed
            D.check_z_slice_size()
            
            % Update graph (will be needed by both labels and data display)
            prevsz = D.graph.z_slicesz;
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
            if ~(strcmp(flag,'chgdata') || (strcmp(flag, 'chg') && ~any(chg_dim == [D.active_dim.x, D.active_dim.y])))
                D.graph.set_ticks()
            end
            
            % Update clipping
            chg_clip = strcmp(flag, 'global') || strcmp(D.clipping.span, 'curview');
            if chg_clip, auto_clip(D, false), end
            
            % Update display
            if fn_ismemberstr(flag, {'global', 'chg_dim'})
                % Reset display
                update_display(D, 'global')
            elseif strcmp(flag, 'chgdata')
                % No change in size, all data need to be redisplayed
                update_display(D, 'chgdata')
            else
                % Smart display update
                if ismember(chg_dim_id, D.internal_dim_id)
                    % changes are within elements (the grid arrangement
                    % remains the same)
                    if fn_ismemberstr(flag, {'perm', 'chg'}) ...
                            || (strcmp(flag, 'all') && D.z_slice.header(chg_dim).n == prevsz(chg_dim))
                        flag = 'chgdata'; % no change in size
                    else
                        flag = 'chgdata&blocksize';
                    end
                    update_display(D, flag)
                elseif ~chg_clip
                    % the grid arrangement changes
                    switch flag
                        case 'chg'  % check this case first because it is this one that occurs when going fast through a list
                            update_display(D, 'chg', chg_dim, e.ind)
                        case {'new', 'remove', 'perm'}
                            update_display(D, flag, chg_dim, e.ind)
                        case 'all'
                            n_cur = size(D.h_transform, chg_dim);
                            n = D.z_slice.sz(chg_dim);
                            if n == n_cur
                                update_display(D, 'chgdata')
                            elseif n > n_cur
                                update_display(D, 'new', chg_dim, n_cur+1:n)
                                update_display(D, 'chg', chg_dim, 1:n_cur)
                            else
                                update_display(D, 'remove', chg_dim, n+1:n_cur)
                                update_display(D, 'chg', chg_dim, 1:n)
                            end
                        case 'chg&new'
                            update_display(D, 'chg', chg_dim, e.ind{1})
                            update_display(D, 'new', chg_dim, e.ind{2})
                        case 'chg&rm'
                            update_display(D, 'chg', chg_dim, e.ind{1})
                            update_display(D, 'remove', chg_dim, e.ind{2})
                        otherwise
                            error('flag ''%s'' is not handled', 'flag')
                    end
                elseif chg_clip
                    % all grid elements need to be updated
                    n_cur = size(D.h_transform, chg_dim);
                    n = D.z_slice.sz(chg_dim);
                    if n==n_cur
                        update_display(D, 'chgdata')
                    elseif n>n_cur
                        update_display(D, 'chgdata')
                        update_display(D, 'new', chg_dim, n_cur+1:n)
                    else
                        update_display(D, 'remove', chg_dim, n+1:n_cur)
                        update_display(D, 'chgdata')
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
            elseif ~strcmp(flag, 'chgdata') && isequal(chg_dim, D.color_dim)
                display_color_legend(D)
            end
        end
        function zoom_change(D, e)
            % update graph positions: if data has changed in size,
            % positioning will be updated upon notification of data change;
            % however if data has not changed in size, positioning needs to
            % be updated here
            if ~e.chg_n_out
                c = disable_listener(D.listeners.ax_siz); %#ok<NASGU> % prevent display update following automatic change of axis position
                D.check_z_slice_size % is z_slice too large for being displayed
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
                D.check_z_slice_size() % is z_slice too large for being displayed
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
