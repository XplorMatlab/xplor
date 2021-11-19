classdef DisplayNavigation < xplr.GraphNode
% The displaynavigation class handles all the callbacks of events (mouse
% click, scroll, etc.) that occur in the graph. 
% These can control the following elements:
% - point selection left button click to move the cross in all visible
%                   dimensions
% - ROI selections  drag with right button to make new selections in
%                   the N.selection_dim dimension(s); use also the Selection
%                   menu on the top and selection context menus that appear
%                   when right-clicking on a displayed selection
% - zoom            zoom can be controlled by different actions, the
%                   affected dimensions are detected automatically:
%                   - with the scroll wheel
%                   - drag the mouse left button and draw rectangle region
%                     where to zoom in
%                   - double-click with the mouse left button to reset zoom
%                   - select the label of the intended dimension to make it
%                     active and make a zooming slider appear, then drag
%                     the edges of the slider or double-click on it to
%                     reset zoom

    properties (SetAccess='private')
        D                                       % parent xplr.viewdisplay
        ha = gobjects
        hf = gobjects
        graph
        cross_center
        cross_data_value                          % data value at cross position
        cross = gobjects                        % display cross selector
        sliders = struct('x', [], 'y', []);        % slider objects
        zoom_filters = struct('x', [], 'y', []);    % connected zoom filters
        point_filters = {};
    end
    properties
        selection_filter         % filter being modified by the selections
        selection_display        % displays of SelectionND: for each selection, 2 lines and a text
        selection_save_file       % name of file for saving current selection
        selection_context       
    end
    properties (SetObservable, AbortSet=true)
        show_cross = true;
        cross_color = [0, 0  0];
        cross_alpha = .5;
        selection_dim_id          % dimensions to which these selections apply, identified by its id
        selection_2d_shape = 'ellipse'; % 'poly', 'free', 'rect', 'ellipse', 'ring', 'segment', 'openpoly', 'freeline'
        selection_round_1d_measure = true;
        selection_show = 'shape+name';
        selection_do_edit = false
        selection_prompt_name = false
        selection_at_most_one = false
    end
    properties (Dependent)
        selection_dim            % dimension(s) to which selections apply, identified by its(their) number
    end
    properties (Dependent, SetAccess='private')
        selection               % list of SelectionND object from the filter
    end

    % Constructor
    methods
        function N = DisplayNavigation(D)
            % parent xplr.viewdisplay object and other external objects
            N.D = D;
            N.ha = D.ha;
            N.hf = D.V.hf;
            N.graph = D.graph;
            
            % buttons
            init_buttons(N)
            
            % cross
            N.display_cross()
            
            % sliders
            init_sliders(N)
            
            % data value display
            init_value_display(N)

            % connect sliders to the active dimensions of the display
            % (note that this is in fact redundant with call in
            % viewdisplay.slicechange when viewdisplay object is created)
            connect_zoom_filter(N)

            % axes_click actions
            set(D.ha, 'buttondownfcn', @(u,e)axes_click(N))

            % scroll wheel zooming
            brick.scrollwheelregister(D.ha, @(n)N.wheel_scroll(n))
            
            % selection menu
            uimenu(N.hf, 'Label', 'Selection', 'callback', @(m, e)N.selection_menu(m))
        end
        function init_buttons(N)
        % function init_buttons(N)
        % 3 buttons that control clipping
            
            % first button to adjust clipping with axes_click movements:
            % display image on it indicating how image luminance and
            % contrast change upon axes_click movements
            [ii, jj] = meshgrid(-13:0,13:-1:0);
            x=(0-ii)./(jj-ii) - .5;
            x(end)=0;
            u = uicontrol('parent', N.D.hp, ...
                'enable', 'inactive', 'cdata', brick.clip(sin(pi*x), [-1, 1], 'gray'), ...
                'buttondownfcn', @(u,e)move_clip(N));
            brick.controlpositions(u, N.ha,[1, 1, 0, 0],[-1, -16, 16, 16])

            % two next buttons control extent of clipping
            u = uicontrol('parent', N.D.hp, ...
                'string', '+', 'fontsize', 8, ...
                'callback', @(u,e)clip_change_range(N, '+'));
            brick.controlpositions(u, N.ha, [1, 1, 0, 0], [-1, -32, 16, 16])
            u = uicontrol('parent', N.D.hp, ...
                'string', '-', 'fontsize', 8, ...
                'callback', @(u,e)clip_change_range(N, '-'));
            brick.controlpositions(u, N.ha, [1, 1, 0, 0], [-1, -48, 16, 16])
        end
        function init_sliders(N)
            N.sliders.x = brick.slider('parent', N.D.hp, 'mode', 'area', ...
                'layout', 'right', 'callback', @(u, evnt)chg_zoom(N, 'x', u));
            N.sliders.y = brick.slider('parent', N.D.hp, 'mode', 'area', ...
                'layout', 'down', 'callback', @(u, evnt)chg_zoom(N, 'y', u));
            pcol = get(N.D.hp, 'backgroundcolor');
            set([N.sliders.x, N.sliders.y], 'visible', 'off', 'scrollwheel', 'on', 'value', [0, 1], ...
                'backgroundcolor', pcol*.75, 'slidercolor', pcol*.95)
            brick.controlpositions(N.sliders.x, N.ha, [0, 1, 1, 0], [0, 0, 0, 12]);
            brick.controlpositions(N.sliders.y, N.ha, [1, 0, 0, 1], [0, 0, 12, -48]);
        end
        function init_value_display(N)
            N.cross_data_value = uicontrol('Parent', N.D.hp, 'style', 'text', 'enable', 'inactive', ...
                'fontsize', 8, 'horizontalalignment', 'right');
            brick.controlpositions(N.cross_data_value, N.D.hp, [1, 0], [-100, 10, 90, 15])
        end
    end
    
    % Get/Set dependent
    methods
        function d = get.selection_dim(N)
            d = N.D.slice.dimension_number(N.selection_dim_id);
        end
        function set.selection_dim(N, dim)
            N.selection_dim_id = dim; % dim will be concerted to dim_id in N.set.selection_dim_id anyway
        end
    end
    
    % Clip control
    methods
        function move_clip(N)
            switch get(N.hf, 'selectiontype')
                case 'normal'       % change clip
                    clip0 = brick.row(N.D.current_clip);
                    if isempty(clip0)
                        % cross outside of current zoom
                        return
                    end
                    e0 = diff(clip0);
                    clip_center = N.D.clipping.center;
                    if ~isempty(clip_center), clip0 = clip_center + [-.5, .5]*e0; end
                    p0 = get(N.hf, 'currentpoint');
                    ht = uicontrol('style', 'text', 'position', [2, 2, 200, 17], 'parent', N.hf);
                    set(ht, 'string', sprintf('min: %.3f,  max: %.3f', clip0(1), clip0(2)))
                    c = onCleanup(@()delete(ht)); 
                    % change clip
                    brick.buttonmotion(@move_clipsub, N.hf, 'pointer', 'cross')
                case 'open'         % use default clipping
                    N.D.auto_clip(~N.D.clipping.buttons_only_for_current_cells)
            end
            function move_clipsub
                % 'naive' new clip
                p = get(N.hf, 'currentpoint');
                dp = p - p0;
                if ~isempty(clip_center), dp = [-1, 1]*(dp(2)-dp(1))/2; end
                % (tolerance to detect motion)
                tol = 2;
                dp = dp - (sign(dp) .* min(tol, abs(dp)));
                FACT = 1/100;
                clip = clip0 + dp*(e0*FACT);
                % it might be that we have diff(clip)<=0 here! apply some
                % transformation to solve that
                e = diff(clip);
                thr = e0/10;
                if e<thr
                    %e = thr*exp(e/thr-1); % goes from thr for e=thr to 0 for e=-Inf
                    e = thr^2 / (2*thr-e); % goes from thr for e=thr to 0 for e=-Inf
                    clip = mean(clip) + [-.5, .5]*e;
                    if clip(1) < clip0(1)
                        clip = clip0(1) + [0 1] * e;
                    elseif clip(2) > clip0(2)
                        clip = clip0(2) + [-1 0] * e;
                    end
                end     
                % update display
                set(ht, 'string', sprintf('min: %.3f,  max: %.3f', clip(1), clip(2)))
                N.D.set_clip(clip, false)
            end
        end
        function clip_change_range(N, flag)
            % current clip extent
            if N.D.clipping.buttons_only_for_current_cells
                clip = N.D.current_clip;
            else
                clip = N.D.grid_clip; % ND array
            end
            m = mean(clip);
            e = diff(clip);

            % round it to a nice value
            e10 = 10.^floor(log10(e));
            e = e ./ e10;
            vals = [.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 15];
            f = sum(bsxfun(@gt, e(:)*1.1 ,vals), 2);
            
            % update as specified
            f = f + brick.switch_case(flag, '+', -1, '-', 1);
            e = e10 .* reshape(vals(f), size(e));
            clip = brick.add(m, brick.mult([-.5; .5], e));
            
            % set clip
            if N.D.clipping.buttons_only_for_current_cells
                N.D.set_clip(clip)
            else
                N.D.set_grid_clip(clip)
            end
        end
    end
    
    % general mouse actions inside the graph
    methods
        function axes_click(N, flag)
            % function axes_click(N)
            % function axes_click(N, 'point_only')
            %---
            % Options:
            % - 'point_only' set this flag when we have already clicked
            %               and released the axes_click button (for example
            %               because we clicked on the cross): in this case
            %               do not start drag zooming or ROI selection
            
            if nargin<2, flag = ''; end
            point_only = strcmp(flag, 'point_only');
            point =  get(N.D.ha, 'CurrentPoint');
            point = point(1, [1, 2])';
            switch get(N.hf, 'SelectionType')
                case 'normal'
                    % zoom in or select point
                    if point_only
                        do_zoom = false;
                    else
                        rect = brick.mouse(N.ha, 'rectaxp-');
                        do_zoom = any(any(abs(diff(rect, 1, 2))>1e-2));
                    end
                    if do_zoom
                        % apply some rules to guess in which dimension the
                        % user intended to zoom, and which value he wanted
                        % to zoom in
                        [zoom_dim, zoom] = N.zoom_in_rect(rect);
                        N.D.zoom_slicer.set_zoom(zoom_dim, zoom)
                    else
                        N.manual_click_move_cross(point);
                    end
                case 'open'
                    % zoom reset
                    % (automatic determination of the intended dimensions
                    % for changing zoom)
                    zoom_dim = N.pointed_zoom_dimension(point, true);
                    % (zoom reset)
                    zoom = repmat(':', 1, length(zoom_dim));
                    N.D.zoom_slicer.set_zoom(zoom_dim, zoom)
                case 'alt'
                    % on right click: create a new selection in N.selection
                    % depending on the parameter selectionshape
                    if ~isempty(N.selection_dim_id)
                        sel_slice = N.selection_mouse();
                        if isempty(sel_slice), return, end

                        % prompt for selection name
                        options = {};
                        if N.selection_prompt_name
                            name = inputdlg('Selection name', 'xplor');
                            if ~isempty(name), options = {'Name', name{1}}; end
                        end
                        
                        % update filter
                        if ~N.selection_at_most_one || isempty(N.selection)
                            N.selection_filter.update_selection('new', sel_slice, options{:})
                        elseif isscalar(N.selection)
                            N.selection_filter.update_selection('chg', 1, sel_slice, options{:})
                        else
                            N.selection_filter.update_selection('chg&rm', {1, 2:length(N.selection)}, sel_slice, options{:})
                        end
                    end
            end
        end
    end
    
    % Cross point selection
    methods
        function connect_point_filter(N, dim, key)
            if nargin < 3, key = 1; end

            % Disconnect current point filters
            if nargin < 2
                disconnect_point_filter(N)
                dim = 1:N.D.slice.nd;
            else
                disconnect_point_filter(N, dim)
            end
            
            % Connect new point filters: loop on dimensions
            for d = dim
                link_key = key;
                head = N.D.slice.header(d);
                % no interest in creating and controling a filter for a
                % dimension with only 1 value
                if head.n == 1
                    N.point_filters{d} = [];
                    continue
                end
                % get filter from bank or create one for the header in this
                % dimension
                P = xplr.Bank.get_point_filter(link_key, head, N); % FilterAndPoint filter
                N.point_filters{d} = P;
                % listen to the point filter event
                N.add_listener(P, 'changed_operation', @(u,e)moved_point(N, e))
            end
        end
        function disconnect_point_filter(N, dim)
            if nargin < 2
                dim = 1:length(N.point_filters);
            end
            for d = dim
                P = N.point_filters{d};
                if isempty(P), continue, end
                xplr.Bank.unregister_filter(P, N)
                N.disconnect(P)  % this is the same as F.disconnect(N)!
            end
        end
        function moved_point(N, e)
            if ~strcmp(e.type, 'point'), return, end
            ijk = get_point_index_position(N);
            N.cross_center = N.graph.slice_to_graph(ijk);
            update_cross_visibility(N);
        end
        function ijk = get_point_index_position(N,varargin)
            % function ijk = getPointIndexPosition(N[,subdimflag][,'global|round'])
            %---
            % Position in slice of point selected by the point filters.
            % Argument subdimflag serves to subselect dimensions of
            % interest and is any of the following: 'all'
            % [default],'internal','external','overlap','overlap+external'.
            % Values for non-selected dimensions will be replaced by ones.
            % If 'global' flag is set, returns global index instead of
            % per-dimensions coordinates. Warning: this global index is to
            % be considered inside an array with only the selected
            % dimensions.
            %
            % Example:
            % slice size is 10*5*4*3, display mode 'time courses',
            % dimension positions {'x', 'mergeddata', 'y', 'x'}, current
            % point [8.1; 5; 3; 2]
            % getPointIndexPosition(N)  returns [8.1; 5; 3; 2]
            % getPointIndexPosition(N,'internal')   returns 8.1
            % getPointIndexPosition(N,'external')   returns [3; 2]
            % getPointIndexPosition(N,'overlap')    returns 5
            % getPointIndexPosition(N,'overlap+external')   returns [5; 3; 2]
            % getPointIndexPosition(N,'clip')       returns [5; 3; 2]
            % getPointIndexPosition(N,'round')      returns [8; 5; 3; 2]
            % getPointIndexPosition(N,'global')     returns 348
            % getPointIndexPosition(N,'internal','global')  returns 8
            % getPointIndexPosition(N,'external','global')  returns 7

            nd = N.D.slice.nd;
            ijk = ones(nd, 1);

            % Which dimensions to consider (note that in every case they
            % will be sorted)
            dim = 1:nd; do_global = false; do_round = false;
            for i = 1:length(varargin)
                switch varargin{i}
                    case 'all'
                        dim = 1:nd;
                    case 'internal'
                        dim = N.D.internal_dim;
                    case 'external'
                        dim = N.D.external_dim;
                    case 'overlap'
                        dim = N.D.overlap_dim;
                    case 'overlap+external'
                        dim = unique([N.D.overlap_dim, N.D.external_dim]);
                    case 'clip'
                        dim = unique(N.D.clip_dim);
                    case 'global'
                        do_round = true;
                        do_global = true;
                    case 'round'
                        do_round = true;
                    otherwise
                        error('unknown flag')
                end
            end
            for d = dim
                P = N.point_filters{d};
                if isempty(P), continue, end
                ijk(d) = P.index_exact;
            end

            % Global index?
            if do_round
                ijk = brick.coerce(round(ijk), 1, N.D.slice.sz');
            end
            if do_global
                ijk = brick.indices(N.D.slice.sz(dim), ijk(dim), 'i2g');
            end
        end
        function reposition_cross(N)
            ijk = get_point_index_position(N);
            N.cross_center = N.graph.slice_to_graph(ijk);
            update_cross_visibility(N);
        end
        function display_cross(N)
           
            % cross
            N.cross(1) = line('Parent', N.D.ha, 'ydata', [-.5, .5]);
            N.cross(2) = line('Parent', N.D.ha, 'xdata', [-.5, .5]);
            N.cross(3) = line('Parent', N.D.ha, 'xdata', [0, 0], 'ydata', [0, 0]); % a single point
            set(N.cross, 'Color', [N.cross_color, N.cross_alpha]) % cross is semi-transparent!
            
            % position
            N.cross_center = [0, 0];
            
            % callbacks
            for i=1:3
                set(N.cross(i),'buttondownfcn', @(u,e)manual_move_cross(N, i))
            end
        end
        function set.cross_center(N, cross_center)

            % re-display cross if it was deleted (happens upon
            % D.resetDisplay)
            if ~all(isvalid(N.cross))
                brick.delete_valid(N.cross)
                N.display_cross()
            end
            
            % set property
            N.cross_center = cross_center;

            % move the cross
            set(N.cross(1), 'XData', cross_center([1, 1]))
            set(N.cross(2), 'YData', cross_center([2, 2]))
            set(N.cross(3), 'XData', cross_center([1, 1]), 'YData', cross_center([2, 2]))

        end
        function manual_move_cross(N, il)
            if ~ismember(get(N.hf, 'selectiontype'), {'normal', 'open'})
                % not a left click: execute callback for axes
                axes_click(N)
                return
            end
            
            % prepare a time to pan zoom while moving the cross!
            do_drag_zoom = (il == 1) && isequal(N.D.layout.x, 1) && isequal(N.D.active_dim.x, 1);
            if do_drag_zoom
                Z = N.D.zoom_filters(1);
                do_drag_zoom = ~strcmp(Z.zoom, ':');
            end
            if do_drag_zoom
                drag_timer = timer('timerfcn', @drag_zoom, 'ExecutionMode', 'fixedSpacing', 'period', .01);
                zoom_width = diff(Z.zoom);
                zoom_min = .5;
                zoom_max = .5 + Z.header_in.n;
                drag_zone = .15; % width of special zone where zoom dragging occurs
                drag_speed = 0;
            end
            point = [];
            
            pointer = brick.switch_case(il, 1, 'left', 2, 'top', 3, 'crosshair');
            anymove = brick.buttonmotion(@move_cross_sub, N.hf, 'moved?', 'pointer', pointer);
            if do_drag_zoom
                stop(drag_timer)
            end
            if ~anymove
                % execute callback for axes
                axes_click(N, 'point_only')
                return
            end
            
            function move_cross_sub
                %anymove = true;
                p = get(N.D.ha, 'currentpoint');
                p = p(1, 1:2);
                switch il
                    case 1
                        point = [p(1), N.cross_center(2)];
                    case 2
                        point = [N.cross_center(1), p(2)];
                    case 3
                        point = p;
                end
                % move the cross
                manual_click_move_cross(N, point)
                % pan view if we are close to edge!
                if do_drag_zoom
                    if abs(p(1)) > .5 - drag_zone
                        drag_speed = sign(p(1)) * ((abs(p(1)) - (.5-drag_zone))/drag_zone)/10;
                        if ~brick.boolean(drag_timer.Running)
                            start(drag_timer)
                        end
                    else
                        stop(drag_timer)
                    end
                end
            end
            
            function drag_zoom(~, ~)
                zoom = Z.zoom + drag_speed * zoom_width;
                if zoom(1) < zoom_min
                    zoom = zoom_min + [0, zoom_width];
                elseif zoom(2) > zoom_max
                    zoom = zoom_max - [zoom_width, 0];
                end
                Z.set_zoom(zoom)
                manual_click_move_cross(N, point)
            end
            
        end        
        function manual_click_move_cross(N, point)
            % move the cross to the selected point
            ijk = brick.row(N.graph.graph_to_slice(point, 'invertible', true));
            
            % clip to data edges
            ijk = brick.coerce(ijk, .5001, N.D.slice.sz + .4999);
            
            % round indices values in dimensions with categorical headers
            categorical = [N.D.slice.header.categorical];
            ijk(categorical) = round(ijk(categorical));
            
            % update the point filters (only for dimensions where the point
            % remains within the slice)
            for d = 1:N.D.nd
                P = N.point_filters{d};
                if ~isempty(P)
                    P.index_exact = ijk(d);
                end
            end
        end        
        function remove_cross(N)
            delete(N.cross)
        end      
        function update_cross_visibility(N)
            % do not show cross?
            if ~N.show_cross
                set(N.cross, 'Visible', 'off')
                set(N.cross_data_value, 'Visible', 'off')
                return
            end           
                        
            % dims that are out of display
            ijk = get_point_index_position(N);
            dim_out_of_display = N.is_index_out_of_display(ijk, true);
           
            % Hide the vertical bar if all dimensions on x are singletons or if
            % cross_center is out of display on any dimension on x
            layout = N.D.layout;
            x_dim = [layout.x, layout.xy];
            x_singleton = isempty(x_dim);
            x_is_out_of_display = any(dim_out_of_display(x_dim));
            N.cross(1).Visible = brick.onoff(~x_singleton && ~x_is_out_of_display);
            
            % Same things for horizontal bar
            y_dim = [layout.y, layout.xy];
            y_singleton = isempty(y_dim);
            y_is_out_of_display = any(dim_out_of_display(y_dim));
            N.cross(2).Visible = brick.onoff(~(y_singleton|y_is_out_of_display));
            
            % Cross center
            update_cross_center_visibility(N);

            % Cross Value
            set(N.cross_data_value,'Visible','on')
            update_value_display(N);
        end
        function update_cross_center_visibility(N)
            %  if one of the dimension of the cross is hidden, hide the
            % cross center as well
            N.cross(3).Visible = brick.onoff(brick.boolean(N.cross(1).Visible) && brick.boolean(N.cross(2).Visible));
        end
        function update_value_display(N)
            idx = get_point_index_position(N, 'global');
            value = N.D.slice.data(idx);

            % Test to display the value as "val(d1,d2,d3,...)=value"
            %set(N.cross_data_value,'String',['val(' num2str(ijk(1),'%.3g') ',' num2str(ijk(2),'%.3g') ')=' ...
            %            num2str(value,'%.3g')])         
            set(N.cross_data_value, 'String', ['Value: ', num2str(value, '%.3g')])
        end
        % cross color, transparency, and global visibility
        function set.show_cross(N, value)
            N.show_cross = value;
            N.update_cross_visibility()
        end
        function set.cross_color(N,color)
            color = brick.colorbyname(color);
            if isempty(color), error 'wrong color', end
            N.cross_color = color;
            set(N.cross, 'color', [N.cross_color, N.cross_alpha])
        end
        function set.cross_alpha(N, alpha)
            N.cross_alpha = alpha;
            set(N.cross, 'color', [N.cross_color, N.cross_alpha])
        end
    end
    
    % ROI selection
    methods
        function selection_menu(N, m)
            % Selection menu is populated when being activated
            delete(get(m, 'children'))
            
            header = N.D.slice.header;
            sel_dim = N.selection_dim;
            sel_labels = {header(sel_dim).label};
            org = N.D.layout;

            % Set selection dimension
            % (no active control? -> propose selection in the 'internal' dimensions)
            prop_dim = [];
            if isempty(N.selection_dim_id)
                if ~isempty(org.xy) && strcmp(org.xy_mode, 'map')
                    % map -> 1D selection in map
                    prop_dim = org.xy;
                elseif strcmp(N.D.display_mode, 'time courses')
                    % time courses -> 1D vertical selection
                    if ~isempty(org.x), prop_dim = org.x(1); end
                else
                    if ~isempty(org.x) && ~isempty(org.y)
                        % image -> 2D selection
                        prop_dim = [org.x(1), org.y(1)];
                    elseif ~isempty(org.x)
                        % image but no Y dimension -> 1D horizontal selection
                        prop_dim = org.x(1);
                    elseif ~isempty(org.y)
                        % image but no x dimension -> 1D vertical selection
                        prop_dim = org.y(1);
                    end
                end
            end
            if ~isempty(prop_dim)
                prop_labels = {header(prop_dim).label};
                if isscalar(prop_dim)
                    item_label = ['Control selection in dimension ', prop_labels{1}];
                else
                    item_label = ['Control selection in dimensions ', brick.strcat(prop_labels,',')];
                end
                prop_dim_id = [header(prop_dim).dim_id];
                uimenu(m, 'label', item_label, ...
                    'callback', @(u,e)set(N, 'selection_dim_id', prop_dim_id))
            end
            % (sub-menu with other possible dimensions)
            switch length(N.selection_dim_id)
                case 0
                    info = 'Control selection in dimension...';
                case 1
                    info = ['Control selection in dimension: ', sel_labels{1}];
                case 2
                    info = ['Control selection in dimensions: ', brick.strcat(sel_labels,',')];
                otherwise
                    error 'programming: selection control in more than 2 dimensions'
            end
            % (submenu for other possibilities)
            m2 = uimenu(m, 'label', info);
            % (1D: dimension location must be x, y or xy)
            dim_ok = sort([org.x, org.y, org.xy]);
            brick.propcontrol(N, 'selection_dim_id', ...
                {'menugroup', {header(dim_ok).dim_id}, {header(dim_ok).label}}, ...
                {'parent', m2});
            % (2D: dimension locations must be respectively x and y)
            available = cell(2, length(org.x), length(org.y));
            if ~isempty(available)
                % selections with first dim on x-axis, second dim on y-axis
                for i = 1:length(org.x)
                    for j = 1:length(org.y)
                        d = [org.x(i), org.y(j)];
                        available{1, i, j} = [header(d).dim_id];
                        available{2, i, j} = brick.strcat({header(d).label}, ',');
                    end
                end
                brick.propcontrol(N, 'selection_dim_id', ...
                    {'menugroup', available(1, :) available(2, :)}, ...
                    {'parent', m2});
                % selections with first dim on y-axis, second dim on x-axis
                for i = 1:length(org.x)
                    for j = 1:length(org.y)
                        d = [org.y(j), org.x(i)];
                        available{1, i, j} = [header(d).dim_id];
                        available{2, i, j} = brick.strcat({header(d).label}, ',');
                    end
                end
                brick.propcontrol(N, 'selection_dim_id', ...
                    {'menugroup', available(1, :), available(2, :)}, ...
                    {'parent', m2});
            end

            %             uimenu(m2,'label','ND selection...','separator','on', ...
            %                 'callback',@(u,e)set(N,'selection_dim_id','prompt'))
            % (stop)
            if ~isempty(N.selection_dim_id)
                uimenu(m, 'label', 'Stop selection control', ...
                    'callback', @(u,e)set(N, 'selection_dim_id', []))
            end

            % Options below should not be displayed if there is no active
            % selection control
            if isempty(N.selection_dim_id), return, end

            % Selection options
            uimenu(m, 'label', 'Clear selections', 'separator', 'on', ...
                'callback', @(u,e)N.selection_filter.update_selection('reset'))
            if length(N.selection_dim_id) == 2
                brick.propcontrol(N, 'selection_2d_shape', ...
                    {'menuval', {'poly', 'free', 'rect', 'ellipse', 'ring', 'line', 'openpoly', 'freeline'}}, ...
                    'parent', m, 'label', 'Shape');
            end
            if length(N.selection_dim_id) == 1 && N.D.slice.header(sel_dim).is_measure
                brick.propcontrol(N, 'selection_round_1d_measure', 'menu', ...
                    {'parent', m, 'label', 'Round 1D selections to data indices', 'separator', 'on'});
            end
            brick.propcontrol(N, 'selection_prompt_name', 'menu', ...
                {'parent', m, 'label', 'Prompt for name of new selections'});
            brick.propcontrol(N, 'selection_at_most_one', 'menu', ...
                {'parent', m, 'label', 'At most one selection'})
            
            brick.propcontrol(N, 'selection_show', ...
                {'menuval', {'shape+name', 'shape', 'name', 'center'}}, ...
                {'parent', m, 'label', 'Display mode', 'separator', 'on'});
            % selection edit mode not ready yet!!! (need first to convert
            % selection from slice to graph)
            %             brick.propcontrol(N,'selection_do_edit','menu', ...
            %                 {'parent',m,'label','Selections modifyable'});

            % Load/save selections
            uimenu(m, 'label', 'Load...', 'separator', 'on', ...
                'callback', @(u,e)N.selection_load())
            uimenu(m, 'label', 'Save', 'enable', brick.onoff(~isempty(N.selection_save_file)), ...
                'callback', @(u,e)N.selection_save(N.selection_save_file))
            uimenu(m, 'label', 'Save as...', ...
                'callback', @(u,e)N.selection_save())
        end
        function sel_slice = selection_mouse(N, message)
            % function sel_slice = selection_mouse(N [,message])
            %---
            % 'pointaction' argument is a function handle that sets a
            % specific function to execute when drag action occurs to stay
            % on a point: this is used to raise a context menu when
            % right-clicking on a selection
            %
            % if second argument 'message' is set, this message will be
            % displayed under the mouse cursor, and it is assumed that the
            % first point of the selection was not clicked yet; otherwise
            % it is assumed that the first point of the selection was
            % already clicked
            
            sel_dim = N.selection_dim;
            selnd = length(sel_dim);
            sz = N.D.slice.sz(sel_dim); % size of data in the dimension where selection is made

            % check dimension location
            org_id = N.D.layout_id;
            dim_location = brick.strcat(org_id.dim_locations(sel_dim), ',');
            if ~ismember(dim_location, {'x', 'y', 'xy', 'x,y', 'y,x'})
                disp(['selection in location ''' dim_location ''' not handled'])
                sel_slice = [];
                return
            end
            use_map = strcmp(dim_location, 'xy') && strcmp(org_id.xy_mode, 'map');
            
            % two different behaviors for the brick.mouse function
            if nargin<2
                % we assume the first point of the selection was already
                % clicked
                popt = '-';
                msgopt = {};
            else
                % we assume the first point of the selection was not
                % clicked yet, and a message will be displayed
                popt = '';
                msgopt = {message};
            end

            % define selection
            if selnd == 1 && ~use_map
                % user interaction
                if isscalar(dim_location)
                    poly_ax = brick.mouse(N.ha, [dim_location, 'segment', popt], msgopt{:});
                    switch dim_location
                        case 'x'
                            poly_ax = [poly_ax; 0, 0];
                        case 'y'
                            poly_ax = [0, 0; poly_ax];
                    end
                else
                    poly_ax = brick.mouse(N.ha, ['rectaxp', popt], msgopt{:});
                end

                % convert to slice indices in the selected
                % dimension
                ijk0 = round(N.graph.graph_to_slice(poly_ax(:, 1)));
                poly_slice = N.graph.graph_to_slice(poly_ax, 'sub_dim', sel_dim, 'ijk0', ijk0);
                poly_slice = sort(poly_slice);

                % create selection in slice indices coordinates
                if N.D.slice.header(sel_dim).categorical
                    sel_slice = xplr.SelectionND('indices', round(poly_slice(1)):round(poly_slice(2)), sz);
                elseif diff(poly_slice)==0
                    sel_slice = xplr.SelectionND('point1D', poly_slice(1));
                elseif N.selection_round_1d_measure
                    sel_slice = xplr.SelectionND('line1D', round(poly_slice) + [-.5, .5]);
                else
                    sel_slice = xplr.SelectionND('line1D', poly_slice);
                end
            elseif selnd == 2 || (selnd == 1 && use_map)
                % user interaction
                mouse_sel_mode = brick.switch_case(N.selection_2d_shape, ...
                    'line', 'segment', {'poly', 'openpoly'}, 'polypt', 'freeline', 'free', ...
                    'ellipse', 'ellipse*', 'ring', 'ring*', ...
                    N.selection_2d_shape);
                sel_type = brick.switch_case(N.selection_2d_shape, ...
                    {'poly', 'free'}, 'poly2D', 'rect', 'rect2D', ...
                    'ellipse', 'ellipse2D', 'ring', 'ring2D', ...
                    {'line', 'openpoly', 'freeline'}, 'openpoly2D');
                poly_ax = brick.mouse(N.D.ha, [mouse_sel_mode popt], msgopt{:});

                % create selection in graph coordinates
                sel_ax = xplr.SelectionND(sel_type, poly_ax);
                % if selection is too small, convert it to a single point
                sel_ax.check_point(.005)

                % convert to slice coordinates
                sel_slice = N.graph.selection_to_slice(sel_dim, sel_ax);
            end
        end
        function sel = get.selection(N)
            F = N.selection_filter;
            if isempty(F)
                sel = [];
            else
                sel = F.selection;
            end
        end
        function set.selection_dim_id(N, dim_id)
            % function setselection_dim_id(N,dim_id)
            %---
            % Select the dimension(s) for which selections are displayed and
            % made in the display
            
            % special: prompt user for selecting dimensions
            if ischar(dim_id) && strcmp(dim_id, 'prompt')
                header = N.D.slice.header;
                sel_dim = listdlg( ...
                    'PromptString', 'Select up to 2 dimensions', ...
                    'ListString', {header.label});
                dim_id = [header(sel_dim(1:min(2,end))).dim_id];
            end

            % check dimension(s)
            [dim, dim_id] = N.D.slice.dimension_number_and_id(dim_id);
            nd = length(dim_id);
            if ~ismember(nd, [0, 1, 2])
                error 'number of dimension for selection display must be 1 or 2'
            end
            if iscell(dim), error 'some dimension is not present in the slice data', end
            singleton = (N.D.slice.sz(dim) == 1);
            dim_id(singleton) = []; % no selection in singleton dimension(s)
            dim(singleton) = [];
            if isequal(dim_id, N.selection_dim_id)
                return
            end
            
            % set property
            N.selection_dim_id = dim_id;
            
            % disconnect from previous filter
            F = N.selection_filter;
            if ~isempty(F)
                xplr.Bank.unregister_filter(F,N)
                N.disconnect(F)
            end
            
            % no selection?
            if isempty(dim_id)
                N.selection_filter = [];
                N.display_selection()
                return
            end
            
            % find filter to connect to
            header_in = N.D.slice.header(dim);
            F = xplr.Bank.get_filter_filter(1, header_in, N); % xplr.filter object
            if isempty(F)
                % filter needs to be created
                error 'not implemented yet'
            end
            N.selection_filter = F;
            
            % watch changes in filter to update display!
            N.add_listener(F, 'changed_operation', @(u,e)selection_filter_change(N, e));
            
            % update display
            N.display_selection()
            
        end
        function check_selection_filter(N)
            % Check whether current dimension for selections display is
            % still valid, i.e. whether the connected filter still fits a
            % dimension in the new slice
            if isempty(N.selection_dim_id), return, end
            if ~all(ismember(N.selection_dim_id, [N.D.slice.non_singleton_header().dim_id]))
                N.selection_dim_id = [];
            end
        end
        function selection_filter_change(N, e)
            % Update selection display
            if ~strcmp(e.type, 'filter'), return, end
            N.display_selection(e.flag, e.ind)
        end
        function display_selection(N, flag, ind)
            % function display_selection(N)
            % function display_selection(N,'reset')
            % function display_selection(N,'all')
            % function display_selection(N,'referentialchanged')
            % function display_selection(N,'new',ind)
            % function display_selection(N,'chg',ind)
            % function display_selection(N,'chg&new',{indchg indnew})
            % function display_selection(N,'chg&rm',{indchg indrm})
            % function display_selection(N,'remove',ind)
            % function display_selection(N,'perm',perm)
           
            if isempty(N.selection_dim_id)
                brick.delete_valid(N.selection_display)
                N.selection_display = [];
                return
            end
            
            % before all, check that selections can be displayed
            selection_axis = [N.D.layout_id.dim_locations{N.selection_dim}];
            if ~ismember(selection_axis, {'x', 'y', 'xy'})
                disp('selections cannot be displayed')
                brick.delete_valid(N.selection_display)
                N.selection_display = [];
                return
            end
            
            % or display update
            if nargin<2, flag = 'all'; end
            if nargin<3, ind = 1:length(N.selection); end
            switch flag
                case {'all', 'reset', 'new'}
                    % delete current display
                    if brick.ismemberstr(flag, {'all', 'reset'})
                        brick.delete_valid(N.selection_display)
                        N.selection_display = [];
                    end
                    % draw new selections
                    for idx=ind, display_one_sel(N, idx, 'new'); end
                    % keep cross above selections, with the cross center
                    % N.cross(3) at the very top
                    try uistack(N.cross([3, 1, 2]), 'top'), end % can fail when D.resetDisplay is invoked
                case {'add', 'chg', 'affinity'}
                    % might be several indices
                    for k=ind, display_one_sel(N, k, 'chg'); end
                case 'referentialchanged'
                    % it is not the positions of selections that have
                    % changed, but the referential of these positions
                    % relative to the main display axes: we need then to
                    % recompute the positions of each selection display
                    % inside this new referential
                    for k = 1:length(N.selection), display_one_sel(N, k, 'chg'), end
                case 'remove'
                    % delete selected selections
                    brick.delete_valid(N.selection_display(ind))
                    N.selection_display(ind) = [];
                    % index (and therefore displayed name) of some other selections have changed
                    for k = min(ind):length(N.selection), display_one_sel(N, k, 'chg'), end
                case 'perm'
                    perm = ind;
                    N.selection_display = N.selection_display(perm);
                    % index (and therefore displayed name) of some selections have changed
                    for k = find(perm~=1:length(N.selection)), display_one_sel(N, k, 'chg'), end
                case 'chg&new'
                    N.display_selection('chg', ind{1})
                    N.display_selection('new', ind{2})
                case 'chg&rm'
                    N.display_selection('chg', ind{1})
                    N.display_selection('remove', ind{2})
                otherwise
                    error('unknown flag ''%s''', flag)
            end
                
        end
        function display_one_sel(N, k, flag, varargin)
            % function display_one_sel(D,k,'new')       - selection at index k is new
            % function display_one_sel(D,k,'chg')       - selection has changed

            % new selection?
            new_sel = brick.switch_case(flag, 'new', true, 'chg', false);
            
            % position
            % selection in slice coordinates
            sel_ij = N.selection(k);
            name = N.selection_filter.header_out.get_item_names{k};
            % convert selection to displayable polygon
            sel_dim = N.selection_dim;
            selection_axis = [N.D.layout_id.dim_locations{sel_dim}];
            if ~ismember(selection_axis, {'x', 'y', 'xy'})
                error('selection cannot be displayed')
            end
            % selection shape
            [polygon, center] = N.D.graph.selection_mark(sel_dim, sel_ij);
            % additional points for when in selection edit mode
            if N.selection_do_edit
                switch sel_ij.type
                    case {'poly2D', 'mixed', 'point2D', 'line2D'} % TODO: not sure about 'point2D'
                        shape_edit_poly = polygon(:, 1:end-1); % the last point is a repetition of the 1st one
                        shape_edit_info = [];
                    case 'rect2D'
                        shape_edit_poly = polygon(:, 1:4); % the 5th point of polygon is a repetition of the 1st one
                        shape_edit_info = [sel_ij.shapes.points', sel_ij.shapes.vectors'];
                    case {'ellipse2D', 'ring2D'}
                        c = sel_ij.shapes.points;
                        u = sel_ij.shapes.vectors;
                        e = sel_ij.shapes.logic;
                        shape_edit_poly = [c-u, c+u];
                        shape_edit_info = {c, u, e};
                    otherwise
                        error programming
                end
            end

            % name
            name_rotation = brick.switch_case(isscalar(N.selection_dim_id) && length (name)>3, 45, 0);

            % color
            colors = brick.colorset;
            col = colors(mod(k-1, size(colors, 1)) + 1, :);

            % Create / update objects
            if new_sel
                hl = struct('shape', [], 'label', [], 'cross', [], 'handles', []);
                % selection shape
                if strfind(N.selection_show, 'shape')
                    hl.shape = line(polygon(1, :), polygon(2, :), 'Parent', N.D.ha);
                end
                % name
                if strfind(N.selection_show, 'name')
                    hl.label = text(center(1), center(2), name, 'Parent', N.D.ha, ...
                        'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
                        'rotation', name_rotation);
                end
                % center marked with a cross
                if strfind(N.selection_show, 'center')
                    hl.cross = line(center(1), center(2), 'Parent', N.D.ha, ...
                        'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 4);
                end
                % handles to modify selection
                if N.selection_do_edit
                    hl.handles = line(shape_edit_poly(1, :), shape_edit_poly(2, :), 'Parent', N.D.ha, ...
                        'LineStyle', 'none', 'marker', '.', 'UserData', shape_edit_info);
                end
                set(struct_to_array(hl), 'Color', col, ...
                    'ButtonDownFcn', @(u,e)N.selection_clicked(u));
                if isempty(N.selection_display)
                    if k ~= 1, error 'attempting to create N.selection_display at initial index different from 1', end
                    N.selection_display = hl;
                else
                    N.selection_display(k) = hl;
                end            
            else
                hl = N.selection_display(k);
                % update position
                set(hl.label, 'position', center)
                set(hl.shape, 'xdata', polygon(1, :), 'ydata', polygon(2, :))
                set(hl.cross, 'xdata', center(1), 'ydata', center(2))
                if N.selection_do_edit
                    set(hl.handles, 'xdata', shape_edit_poly(1, :), 'ydata', shape_edit_poly(2, :), ...
                        'UserData', shape_edit_info);
                end
                % update name
                set(hl.label, 'string', name, 'rotation', name_rotation)
                set(struct_to_array(hl), 'color', col)
            end
        end
        function set.selection_show(N, value)
            N.selection_show = value;
            % we must see the shape when selection edit is active
            if N.selection_do_edit && ~any(strfind(N.selection_show, 'shape'))
                N.selection_do_edit = false;
            else
                N.display_selection('all')
            end
        end
        function set.selection_do_edit(N, value)
            N.selection_do_edit = value;
            % we must see the shape when selection edit is active
            if N.selection_do_edit && ~any(strfind(N.selection_show,'shape'))
                N.selection_show = 'shape+name';
            else
                N.display_selection('all')
            end
        end
        function set.selection_at_most_one(N, value)
            N.selection_at_most_one = value;
            if value
                nsel = length(N.selection);
                if nsel > 1
                    N.selection_filter.update_selection('remove', 2:nsel)
                end
            end
        end
        function selection_clicked(N, u)
            % click on selection u: 
            % - if selection edit mode is active, some specific actions can
            %   happen
            % - if right-click and no mouse motion, raise a context menu
            % - if none of the two types of actions above are triggered,
            %   execute normal axes callback
            
            % first get index of the clicked selection
            idx = brick.find(@(hl)any(struct_to_array(hl) == u), N.selection_display, 'first');
            
            click_type = get(N.D.V.hf, 'SelectionType');
            if strcmp(click_type, 'alt')
                % perform a regular new selection, but if mouse did not
                % move, show the context menu instead of creating a point
                % selection
                sel_slice = N.selection_mouse();
                if isempty(sel_slice), return, end
                if ~is_point(sel_slice)
                    % prompt for selection name
                    options = {};
                    if N.selection_prompt_name
                        name = inputdlg('Selection name', 'xplor');
                        if ~isempty(name), options = {'Name', name{1}}; end
                    end

                    % update filter
                    if N.selection_at_most_one
                        N.selection_filter.update_selection('chg&rm', {1, 2:length(N.selection)}, sel_slice, options{:})
                    else
                        N.selection_filter.update_selection('new', sel_slice, options{:})
                    end
                else
                    % instead of creating a point selection, raise the
                    % context menu
                    N.selection_context_menu(idx)
                end
            elseif N.selection_do_edit
                flag = '';
                if u == N.selection_display{idx}(1)
                    flag = 'line';
                elseif u == N.selection_display{idx}(end)
                    flag = 'handle';
                end
                if ~isempty(flag)
                    N.selection_do_edit(N, idx, flag)
                else
                    N.axes_click()
                end
            else
                N.axes_click()
            end
        end
        function selection_edit(N, idx, flag)
            % function selection_edit(N,idx,'line|handle|redraw')
            
            switch flag
                case 'redraw'
                    sel_slice = N.selection_mouse('select new shape');
                    N.selection_filter.update_selection('chg', idx, sel_slice)
            end
            
            % code below is not ready
            return
            
            htl = N.selection_display{idx};
            h_shape = hl(1);
            h_handle = hl(end);
            hl = htl(2:end);
            [flag_pt, flag_lin] = brick.flags({'handle', 'line'}, flag);
            p = get(N.ha, 'currentpoint');
            p = p(1, 1:2)';
            polymark = [get(hl(2), 'xdata');
            get(hl(2), 'ydata')];
            sel_dimsnum = N.sel_dims - 'w';
            selection_marks = N.SI.selection.getselset(sel_dimsnum).singleset;
            switch get(N.hf, 'selectiontype')
                case 'normal'               % MOVE POINT
                    shapetype = selection_marks(idx).type;
                    switch shapetype
                        case {'poly2D', 'mixed', 'line2D', 'openpoly2D'}
                            % note that first and last point in polygon are
                            % the same!
                            if flag_pt
                                % closest point
                                dist = sum(brick.add(polymark, -p).^2);
                                [dum, idx] = min(dist);
                                idx = idx(1); %#ok<*ASGLU>
                                if idx==1 && ~brick.ismemberstr(shapetype, {'line2D', 'openpoly2D'})
                                    % need to move both the first and last point (which is a repetition of the first)
                                    idx=[1, size(polymark, s2)];
                                end 
                            else
                                % closest segment (in fact, closest line)
                                a = polymark(:, 1:end-1);
                                b = polymark(:, 2:end);
                                ab = b - a;
                                ab2 = sum(ab.^2);
                                ap = brick.add(p, -a);
                                abap = ab(1, :).*ap(2, :) - ab(2, :).*ap(1, :);
                                dist = abs(abap) ./ ab2;
                                [dum, idx] = min(dist);
                                idx = idx(1);
                                polymark = [a(:, 1:idx), p, b(:, idx:end)];
                                set(hl, 'xdata', polymark(1, :), 'ydata', polymark(2, :))
                                idx = idx + 1;
                            end
                            brick.moveobject(hl, 'point', idx)
                            sel_edit_update_slice(N, idx)
                        case 'rect2D'
                            desc = get_app_data(hl(2), 'description');
                            x = desc(1);
                            y = desc(2);
                            w = desc(3);
                            h = desc(4);
                            x2 = x + w;
                            y2 = y + h;
                            if flag_pt
                                % move corner
                                pol = [x, x2, x2, x; y, y, y2, y2]; % anti-clockwise from (x,y)
                                dist = sum(brick.add(pol, -p).^2);
                                [dum, idx] = min(dist);
                                idx = idx(1);
                            else
                                % move edge
                                dist = abs([p(2)-y, p(1)-x2, p(2)-y2, p(1)-x]);
                                [dum, idx] = min(dist);
                            end
                            col = get(hl(1), 'color');
                            set(hl, 'color', .5*[1, 1, 1])
                            chg_rectangle(N.ha, hl, flag_pt, idx,desc)
                            brick.buttonmotion({@chg_rectangle, N.ha, hl, flag_pt, idx, desc}, N.hf);
                            set(hl, 'color', col)
                            sel_edit_update_slice(N, idx)
                        case {'ellipse2D', 'ring2D'}
                            desc = get_app_data(hl(2), 'description');
                            if flag_pt
                                % closest of two anchor points
                                dist = sum(brick.add(polymark, -p).^2);
                                [dum, idx] = min(dist);
                                idx = idx(1);
                            elseif strcmp(shapetype, 'ellipse2D')
                                % eccentricity
                                idx = 0;
                            else
                                % eccentricity or secondary radius?
                                polygon = [get(hl(1), 'xdata'); get(hl(1), 'ydata')];
                                dist = sum(brick.add(polygon,-p).^2);
                                [dum, idx] = min(dist);
                                if idx < length(polygon)/2
                                    % eccentricity
                                    idx = 0;
                                else
                                    % secondary radius
                                    idx = -1;
                                end
                            end
                            col = get(hl(1), 'color');
                            set(hl, 'color', .5*[1, 1, 1])
                            chg_ellipse(N.ha, hl, idx, desc)
                            brick.buttonmotion({@chg_ellipse, N.ha, hl, idx, desc}, N.hf);
                            set(hl, 'color', col)
                            sel_edit_update_slice(N, idx)
                        otherwise
                            error programming
                    end
                case 'extend'               % MOVE SHAPE
                    if flag_pt
                        dp = brick.moveobject(htl);
                        sel_edit_update_slice(N, idx, dp)
                    elseif flag_lin
                        dp = brick.moveobject([N.seldisp{:}]);
                        sel_edit_update_slice(N, 1:length(N.seldisp), dp) % move all shapes
                    end
                case 'alt'                  % REMOVE
                    if brick.ismemberstr(selection_marks(idx).type, ...
                            {'poly2D', 'mixed'}) && flag_pt
                        % closest point -> remove vertex
                        dist = sum(brick.add(polymark, -p).^2);
                        [dum, idx] = min(dist);
                        idx = idx(1);
                        if idx == 1, idx = [1, size(polymark, 2)]; end
                        polymark(:, idx) = [];
                        sel_ax = SelectionND('poly2D', polymark);
                        sel = AX2IJ(N.SI, sel_ax);
                        update_selection(N.SI, 'change', idx,sel)
                    else
                        % replace the whole shape
                        set(hl, 'visible', 'off')
                        mouse_sel_mode = brick.switch_case(N.shapemode, ...
                            {'poly', 'openpoly'}, 'polypt', 'freeline', 'free', ...
                            N.shapemode);
                        TYPE = brick.switch_case(N.shapemode, {'poly', 'free'}, 'poly2D', 'rect', 'rect2D', ...
                            'ellipse', 'ellipse2D', 'ring', 'ring2D', ...
                            'segment', 'line2D', {'openpoly', 'freeline'}, 'openpoly2D');
                        poly_ax = brick.mouse(N.ha, mouse_sel_mode, 'select new shape');
                        sel_ax = SelectionND(TYPE, poly_ax);
                        sel = AX2IJ(N.SI, sel_ax);
                        update_selection(N.SI, 'change', idx,sel)
                        set(hl, 'visible', 'on')
                    end
            end
        end
        function selection_context_menu(N, idx)
            % init context menu
            if isempty(N.selection_context)
                % create context menu for the first time
                N.selection_context = uicontextmenu(N.D.V.hf);
            end
            m = N.selection_context;
            delete(get(m, 'children'))
            
            % remove selection
            uimenu(m, 'label', 'remove selection', ...
                'callback', @(u,e)N.selection_filter.update_selection('remove', idx))
            
            % replace selection
            uimenu(m, 'label', 're-draw selection', ...
                'callback', @(u,e)N.selection_edit(idx, 'redraw'))
            
            % change selection number
            n_sel = N.selection_filter.n_sel;
            if idx ~= 1
                uimenu(m, 'label', 'make this selection first', ...
                    'callback', @(u,e)N.selection_filter.update_selection('perm', [idx, setdiff(1:n_sel, idx)]))
            end
            if idx ~= n_sel
                uimenu(m, 'label', 'make this selection last', ...
                    'callback', @(u,e)N.selection_filter.update_selection('perm', [setdiff(1:n_sel,idx), idx]))
            end
            
            % make menu visible
            p = get(N.D.V.hf, 'currentpoint');
            p = p(1, 1:2);
            set(m, 'Position', p, 'Visible', 'on')
        end
        function selection_save(N, f_name)
            if nargin<2 || isempty(f_name)
                prompt = 'Select file for saving selections';
                f_name = brick.savefile('*.xpls', prompt, N.selection_save_file);
                if isequal(f_name,0), return, end
            end
            N.selection_filter.save_to_file(f_name);
            N.selection_save_file = f_name;
        end
        function selection_load(N, f_name)
            if nargin < 2
                f_name = brick.getfile('*.xpls', 'Select selections file');
                if isequal(f_name,0), return, end
            end
            try
                N.selection_filter.load_from_file(f_name);
                N.selection_save_file = f_name;
            catch ME
                errordlg(ME.message)
            end
        end
   end

    % Zoom (left mouse, slider and scroll wheel)
    methods
        function zoom_dim = pointed_zoom_dimension(N, point, zooming_out)
            % automatic determination of the intended dimensions for
            % changing zoom: apply to the most internal dimensions where we
            % are not ouf of display
            ok_dim = ~N.is_point_out_of_display(point, true);
            % case of image display: if point is outside image in either of
            % the 2 dimensions, zoom in neither of the 2
            internal_dim = N.D.internal_dim;
            ok_dim(internal_dim) = all(ok_dim(internal_dim));
            % if zooming out, we also want to select dimension(s) where
            % there is an active zoom
            if nargin < 3, zooming_out = false; end
            if zooming_out
                dim_zoom_filters = N.D.zoom_filters;
                ok_dim([dim_zoom_filters.no_zoom]) = false;
                % if no dimension left, consider all dimensions where a
                % zoom is active
                if ~any(ok_dim)
                    ok_dim = ~[dim_zoom_filters.no_zoom];
                end
                % if image display, ok to zoom out in image if ok for any
                % of the 2 dimensions
                ok_dim(internal_dim) = any(ok_dim(internal_dim));
            end
            % most internal x and y dimensions where we are not out of
            % display
            org = N.D.layout;
            zoom_dim = [org.x(find(ok_dim(org.x), 1, 'first')), ...
                org.y(find(ok_dim(org.y), 1, 'first'))];
            % if empty, grid (xy) dimension
            if isempty(zoom_dim)
                zoom_dim = org.xy;
            end
        end
        function [zoom_dim, zoom] = zoom_in_rect(N, rect)
            % determine in which dimension to zoom and zoom values
            %
            % Input:
            % - rect    nd*2 array, first and second point selected by the
            %           mouse
            %
            % Ouput:
            % - zoom_dim    dimension(s) in which to zoom in
            % - zoom        zoom values in this(these) dimension(s)
            
            % since computations rely on the internal representation for
            % dispatching data, they are in a xplr.DisplayGraph method
            [zoom_dim, zslice_zoom, clip_zijk, clip_zoom] = N.D.graph.zoom_in_rect(rect);

            % change in clipping rather than zoom
            if ~isempty(clip_zoom)
                clip = subsref_dim(N.D.grid_clip, 1+N.D.clip_dim, clip_zijk(N.D.clip_dim));
                clip = clip(1) + clip_zoom(:) * diff(clip);
                N.D.set_clip(clip, clip_zijk);
            end
            
            % convert from zslice to slice
            zoom = N.graph.zslice_to_slice(zslice_zoom, false, zoom_dim)'; % 2*nzd
        end
        function chg_zoom(N, f, obj)
            dim = N.D.active_dim.(f);
            if isempty(dim), return, end
            % linked object
            Z = N.D.zoom_filters(dim);
            % prevent unnecessary update of slider display
            c = brick.disable_listener(N.sliders.(f));
            % set value
            if isequal(obj.value, obj.minmax)
                set_zoom(Z,':')
            else
                set_zoom(Z, obj.value)
            end
        end
        function wheel_scroll(N, n_scroll)
%             persistent last_tic last_zoom_dim

            p = get(N.D.ha, 'currentpoint');
            p = p(1, 1:2);
            origin = brick.row(N.graph.graph_to_slice(p)); % current point in data coordinates
            zoom_factor = 1.5^n_scroll;

%             % if we keep zooming from the same point, apply to the same
%             % dimension if it is not the same dimensions that are pointed
%             % any more
%             if ~isempty(last_tic) && toc(last_tic) < .6
%                 zoom_dim = last_zoom_dim;
%             else
%                 % automatic determination of the intended dimensions for
%                 % changing zoom
%                 zoom_dim = N.pointed_zoom_dimension(p);
%             end
%             [last_tic, last_zoom_dim] = deal(tic, zoom_dim);
            zoom_dim = N.D.internal_dim;

            % This commented code had been put to replace the line below,
            % but it seems that the effect is less intuitive. Let's go back
            % to the previous code and see if we get errors or unintuitive
            % behaviors to decide what to do. (TD 12/11/2019)
            %             if n_scroll<0 && ~any(N.D.layout.xy)
            %                 % it does not make sense to zoom-in in a dimensions which
            %                 % does not fill its available space due to aspect ratio
            %                 % constraints
            %                 dim(N.graph.filling(dim)<1) = [];
            %             end
            %             zoom = N.graph.getZoom(dim); %,'effective');
            zoom = N.graph.get_zoom(zoom_dim, 'effective');
            origin_dim = origin(zoom_dim);
            origin_dim = brick.coerce(origin_dim, zoom);
            new_zoom = brick.add(origin_dim, brick.mult(zoom_factor, brick.subtract(zoom, origin_dim)));
            N.D.zoom_slicer.set_zoom(zoom_dim, new_zoom)
        end
    end

    % Update upon changes in active dim and zoom
    methods
        function connect_zoom_filter(N, f)
            % both x and y?
            if nargin < 2
                connect_zoom_filter(N, 'x')
                connect_zoom_filter(N, 'y')
                return
            end
            % slider object and corresponding data dimension
            obj = N.sliders.(f);
            d = N.D.active_dim.(f);
            % disconnect from previous zoom_filters
            z_old = N.zoom_filters.(f);
            if ~isempty(z_old), N.disconnect(z_old), end
            % no active dim?
            if isempty(d)
                set(obj, 'visible', 'off')
                return
            end
            % update slider display to reflect zoom in the specified
            % dimension
            Z = N.D.zoom_filters(d);
            N.zoom_filters.(f) = Z;
            set(obj, 'visible', 'on', 'minmax', [.5, Z.header_in.n + .5], 'value', Z.zoom_value)
            % watch changes in zoom
            function zoom_change(u,e)
                if strcmp(e.type,'zoom')
                    set(obj,'value',Z.zoom_value)
                end
            end
            N.add_listener(Z, 'changed_operation', @zoom_change);
        end
    end
    
    % Tool used by both cross and selection display: is point out of display
    methods
        % return true if point (graph coordinates) is part of the slice
        % data displayed by converting the point to slice coordinates and
        % test if its between minimal and maximal values for all dimensions
        %
        % @param point: 2xn double 
        % @param per_dim: whether to return a single logical value per point
        %                (point is/isn't inside slice, default behavior)
        %                or a vector (for each dimension, if per_dim=true)
        % @return output: 1xn boolean
        function output = is_point_out_of_display(N, point, per_dim)
            if nargin<3, per_dim = false; end
            
            % get slice indices corresponding to the point
            ijk = N.graph.graph_to_slice(point, 'invertible', true)'; % n x ndim
            
            output = is_index_out_of_display(N, ijk, per_dim);
        end
        function output = is_index_out_of_display(N, ijk, per_dim)
            % input
            if nargin<3, per_dim = false; end
            if isvector(ijk) && size(ijk, 2) ~= N.D.slice.nd
                ijk = brick.row(ijk);
            end
            
            % get the min and max slice values of the data displayed
            zoom_slice_values = N.graph.get_zoom('displaylimit'); % 2 x ndim
            
            % set the ouput to zeros (they will be set to one if one of the
            % dimension if it's out of display)
            if per_dim
                output = bsxfun(@lt, ijk, zoom_slice_values(1, :)) ...
                    | bsxfun(@gt, ijk, zoom_slice_values(2, :)); % n x ndim
            else
                output = min(ijk, 1) < zoom_slice_values(1, :) ...
                    | max(ijk, 1) > zoom_slice_values(2, :); % 1 x ndim
            end
        end
    end
    
    
    
end


%--- 
% Matlab's struct_to_array is only provided with some of the toolboxes, so
% let's write our own
function x = struct_to_array(s)

    c = struct2cell(s);
    x = [c{:}];

end
