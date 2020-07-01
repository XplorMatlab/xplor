classdef DisplayLabels < xplr.GraphNode
% displaylabels

    properties (Access='private')
        % parent xplr.viewdisplay object and other external objects
        D
        graph
        % internal
        height     % unit: normalized to axis size
        rot_height  % unit: normalized to axis size
        h
        moving_dim
        moving_clone
        listeners = struct('slice', []);
        prev_org_set_pos     % last organization seen by set_positions
    end
    properties (SetObservable=true)
        do_immediate_display = false;
    end
    properties (Dependent, Access='private')
        ha
    end
        
    methods
        function L = DisplayLabels(D)
            % contructor displaylabels
            
            % parent xplr.viewdisplay object
            L.D = D;
            
            % create labels and position them
            L.graph = D.graph;
            create_labels(L, 'global')
            get_heights(L)
            if ~isempty(D.layout_id), set_positions(L), end
            % watch in axes size (no need to take care of deleting this
            % listener when L is deleted, because there is no reason that L
            % be deleted without D, D.ha, and the listener itself being
            % deleted as well)
            fn_pixelsizelistener(D.ha, @(u,e)update_labels(L, 'axsiz'))
            % note that changes in D lead to display updates through direct
            % method calls rather than notifications
        end
        function ha = get.ha(L)
            ha = L.D.ha;
        end
    end
    
    % Create and position labels
    methods (Access='private')
        function create_labels(L, flag, dims)
            % create labels
            if nargin < 2, flag = 'global'; end
            if nargin < 3 || strcmp(flag, 'global'), dims = 1:L.D.nd; end
            [dims, dim_ids] = L.D.slice.dimension_number_and_id(dims);
            cur_headers = L.D.zslice.header;
            switch flag
                case 'global'
                    delete_valid(L.h)
                    L.h = gobjects(1, L.D.nd);
                otherwise
                    error('invalid flag ''%s''', flag)
            end
            for i = 1:length(dims)
                d = dims(i);
                dim_id = dim_ids(i);
                str = cur_headers(d).label;
                if ~isempty(cur_headers(d).unit), str = [str, ' (', cur_headers(d).unit, ')']; end
                L.h(d) = text('string', ['  ', str, '  '], 'parent', L.ha, ...
                    'margin', 1, ...
                    'backgroundcolor', [1, 1, 1]*.95, 'units', 'normalized', ...
                    'UserData', dim_id, 'buttondownfcn', @(u,e)label_click(L, u), ...
                    'UIContextMenu', uicontextmenu(L.D.V.hf, 'callback', @(m, e)L.D.dimension_context_menu(m, dim_id)));
            end
        end
        function change_label(L, d)
            % change properties that need to when data undergoes a 'chg_dim'
            % change
            
            % label
            head = L.D.zslice.header(d);
            str = head.label;
            if ~isempty(head.unit), str = [str, ' (', head.unit, ')']; end
            set(L.h(d), 'string', ['  ', str, '  '])
            set(L.h(d), 'UserData', head.dim_id)
        end
        function get_heights(L)
            % get heights
            font_size = get(L.D.hp, 'defaultuicontrolfontsize'); % unit: points
            hinch = 1.5*font_size/72;
            hpix = hinch*get(0, 'ScreenPixelsPerInch');
            axsiz = fn_pixelsize(L.ha);
            L.height = hpix/axsiz(2);     % normalized to axis size
            L.rot_height = hpix/axsiz(1);  % normalized to axis size
            L.prev_org_set_pos = []; % next call to set_positions should have doloose set to false
        end
        function set_positions(L, dim)
            % input
            if nargin < 2
                dim = 1:L.D.nd;
            else
                dim = L.D.slice.dimension_number(dim);
            end
            
            % current layout
            org = L.D.layout;

            % locations of dimension: note that some dimensions might not
            % be present in the layout, in particular when it was removed
            % from the layout in labelMove, but the change is not
            % definitive and the slice has not been updated yet
            dim_locations = L.D.layout_id.dim_locations;

            % visible labels and active dimensions
            sz = L.D.slice.sz; % slice size
            ok_dim = (sz > 1) & ~fn_isemptyc(dim_locations);
            isactive = false(1,length(sz));
            isactive([L.D.active_dim.x L.D.active_dim.y]) = true;
            
            % steps in the direction orthogonal to the positioning one
            axis_pos = fn_pixelpos(L.D.ha);
            available_space = axis_pos(1:2)./axis_pos(3:4) - L.height;
            xv_step = min(1.5*L.height, available_space(2)/(1.5+length(org.x)));
            yh_step = min(1.5*L.rot_height, available_space(1)/(1.5+length(org.y)));
            
            % set positions 
            for d = dim
                if ~isgraphics(L.h(d), 'text')
                    % it can happen that labels do not exist yet when
                    % update_labels(L,'axsiz') is invoked
                    continue
                end
                if ~ok_dim(d)
                    % do not display the label of singleton dimension
                    set(L.h(d), 'visible', 'off')
                    if any(d == L.moving_dim)
                        set(L.moving_clone, 'visible', 'off')
                    end
                else
                    f = dim_locations{d};
                    % fixed properties
                    set(L.h(d), 'visible', 'on', ...
                        'rotation', fn_switch(f, {'y', 'yx'}, 90, 0), ...
                        'horizontalalignment', fn_switch(f, {'x', 'y'}, 'center', {'xy', 'merged_data'}, 'left', 'yx', 'right'), ...
                        'verticalalignment', fn_switch(f, 'yx', 'bottom', 'middle'), ...
                        'EdgeColor', fn_switch(isactive(d), 'k', 'none'))
                    % set position
                    switch f
                        case 'x'
                            i = find(org.x == d, 1);
                            x_pos = .5 + L.graph.label_position(d);
                            new_pos = [x_pos, -(1.5+i)*xv_step];
                        case 'y'
                            i = find(org.y == d, 1);
                            y_pos = .5 + L.graph.label_position(d);
                            new_pos = [-(1.5+i)*yh_step, y_pos];
                        case 'merged_data'
                            i = find(org.merged_data == d, 1);
                            new_pos = [-4*L.rot_height, (length(org.merged_data) + .5 - i)*L.height];
                        case 'xy'
                            new_pos = [0, 1 + 1.5*L.height];
                        case 'yx'
                            new_pos = [1+2. 2*L.rot_height, 1];
                    end
                    if any(d == L.moving_dim)
                        set(L.moving_clone, 'position', new_pos, ...
                            'rotation', get(L.h(d), 'rotation'), ...
                            'visible', onoff(true), ...
                            'horizontalalignment', get(L.h(d), 'horizontalalignment'), ...
                            'verticalalignment', get(L.h(d), 'verticalalignment'))
                    else
                        set(L.h(d), 'position', new_pos)
                    end
                end
            end
        end
    end
    
    % Update labels
    methods
        function update_labels(L, flag, chg_dim)
            % input
            if nargin < 2, flag = 'pos'; end
            
            % update
            switch flag
                case 'global'
                    create_labels(L, 'global')
                case 'chg_dim'
                    % some specific properties need to be updated
                    for d = chg_dim, change_label(L, d), end
                case 'axsiz'
                    if isempty(L.D.layout_id) || isempty(L.D.graph.steps)
                        % happens sometimes at init because figure size changes for no clear reason
                        return
                    end 
                    get_heights(L)
                case 'pos'
                case 'active'
                    % only mark labels as active or not
                    nd = L.D.nd;
                    isactive = false(1, nd);
                    isactive([L.D.active_dim.x, L.D.active_dim.y]) = true;
                    for d=1:nd
                        set(L.h(d), 'EdgeColor', fn_switch(isactive(d), 'k', 'none'))
                    end
                    return
                otherwise
                    error('invalid flag ''%s''',flag)
            end
            L.set_positions()
        end
    end
    
    % Label click
    methods
        function label_click(L, obj)
            % dimension number
            dim_id = get(obj, 'UserData');
            % which action
            switch get(L.D.V.hf, 'SelectionType')
                case 'normal'
                    % move the label; if it is not moved, change active dim
                    label_move(L, dim_id)
                case 'open'
                    % position label back in default place
                    L.set_positions(dim_id)
            end
        end
        function label_move(L, dim_id, do_swap, do_initial_move)
            % function label_move(L, dim_id, do_swap, do_initial_move)
            %---
            % Input:
            % - dim_id  identifier of the dimension of the label to move
            % - do_swap if this dimension is in x and is moved to y, allow
            %           swapping with elements that are in y (this results
            %           for example in transposing images)
            % - do_initial_move     label will immediately be positionned
            %           under the mouse pointer
            
            % input
            if nargin<3, do_swap = true; end
            if nargin<4, do_initial_move = false; end

            % prepare for changing organization
            [prev_layout_id, layout_id_d, layout_id] = deal(L.D.layout_id_all); % previous, previous without d, current
            d_layout = layout_id.dim_location(dim_id);
            didx = find(layout_id.(d_layout) == dim_id);
            layout_id_d.(d_layout) = setdiff(layout_id.(d_layout), dim_id, 'stable');
            
            % only one dimension gan go to either xy or yx, and if image
            % display only one dimension can go to merged_data=color_dim
            % -> remember the current dimension in these location for a
            % swap
            xy_dim_id = [layout_id_d.xy, layout_id_d.yx];
            if strcmp(L.D.display_mode, 'image')
                color_dim_id = layout_id_d.merged_data;
            else
                color_dim_id = [];
            end

            % data
            sz = L.D.slice.sz;
            ok_dim = (sz>1);
            dim = L.D.slice.dimension_number(dim_id);
            dim_ids_ok = [L.D.slice.header(ok_dim).dim_id];
            
            % set thresholds: 
            % thresholds will first be computed while ignoring
            % singleton dimensions, then NaNs will be inserted at the
            % locations of these singleton dimensions
            % (x)
            x_layout = layout_id.x;
            x_thr = .5 + L.graph.label_position(x_layout); % first get threshold with d and singleton dimensions still included
            x_thr(x_layout == dim_id) = [];
            x_layout(x_layout == dim_id) = [];        % remove dimension d
            ok_x = ismember(x_layout, dim_ids_ok);
            x_thr(~ok_x) = NaN;             % make it impossible to insert to the right of a singleton or non-present dimension
            if do_swap && strcmp(d_layout, 'y')  % refine intervals in case of swapping
                % not only x insertions, but also x/y swaps are possible
                % for example:
                %      pos =    0   x           x               1
                %   -> thr =       | |       |     |
                nok_x = sum(ok_x);
                pos = [0, x_thr(ok_x), 1]; % add absolute min and max
                x_thr_ok = zeros(2, nok_x);
                for i_thr=1:nok_x
                    dmax = min(pos(i_thr+1) - pos(i_thr), pos(i_thr+2) - pos(i_thr+1))/4; % maximal distance to an x-label to swap with this label
                    x_thr_ok(:, i_thr) = pos(i_thr+1) + [-dmax, dmax];
                end
                x_thr = NaN(2, length(x_layout));
                x_thr(:, ok_x) = x_thr_ok;
            end
            x_thr = [-Inf, row(x_thr)];
            % (y)
            y_layout = layout_id.y;
            y_thr = .5 - L.graph.label_position(y_layout); % first get threshold with d and singleton dimensions still included
            y_thr(y_layout == dim_id) = [];
            y_layout(y_layout == dim_id) = [];        % remove dimension d
            ok_y = ismember(y_layout, dim_ids_ok);
            y_thr(~ok_y) = NaN;             % make it impossible to insert to the right of a singleton or non-present dimension
            if do_swap && strcmp(d_layout, 'x')
                % not only y insertions, but also x/y swaps are possible
                nok_y = sum(ok_y);
                pos = [0, y_thr(ok_y), 1];
                y_thr_ok = zeros(2, nok_y);
                for i_thr=1:nok_y
                    dmax = min(pos(i_thr+1) - pos(i_thr), pos(i_thr+2) - pos(i_thr+1))/4; % maximal distance to an x-label to swap with this label
                    y_thr_ok(:, i_thr) = pos(i_thr+1) + [-dmax dmax];
                end
                y_thr = NaN(2, length(y_layout));
                y_thr(:, ok_y) = y_thr_ok;
            end
            y_thr = [-Inf, row(y_thr)];
                        
            % moving out of the graph to the "filter" area
            view_control = L.D.V.C; % what a nice piece of code isn't it!?
            do_filter = false;

            % label object, make sure it will not be covered by data display
            L.moving_dim = L.D.slice.dimension_number(dim_id);
            obj = L.h(L.moving_dim);
            uistack(obj, 'top')

            % prepare clone
            L.moving_clone = copyobj(obj, L.ha);
            set(L.moving_clone, 'color', [1, 1, 1]*.6, 'edgecolor', 'none', 'BackgroundColor', 'none')
            uistack(L.moving_clone, 'bottom')
            L.set_positions(L.moving_dim)

            % move
            if do_initial_move
                move_label()
            end
            any_change = false;
            moved = fn_buttonmotion(@move_label, L.D.V.hf, 'moved?', 'pointer', 'hand');

            function move_label
                p = get(L.ha,'currentpoint');
                p = p(1, 1:2)';
                p_fig = fn_coordinates(L.ha, 'a2p', p, 'position'); % convert to 'normalized' unit in parent panel
                p = fn_coordinates(L.ha, 'a2c', p, 'position');    % convert to 'normalized' unit in axes
                % move object
                set(obj, 'position', p)
                % special: going outside of graph, in the controls area
                % -> filter dimension, new layout has no more dimension dimID
                if p_fig(1) < 0
                    if ~do_filter
                        do_filter = true;
                        view_control.show_inoperant_filter(dim_id)
                    end
                else
                    if do_filter
                        do_filter = false;
                        view_control.remove_inoperant_filter(dim_id)
                    end
                end
                % immediate display update
                immediate_display = L.do_immediate_display;
                % update organization and object location
                new_layout_id = layout_id_d;
                if p_fig(1) < -2
                    % nothing to do: newlayout_id is already the current
                    % layout without dimension d
                    % however we cannot perform a complete display update
                    % because the slice is not recomputed yet
                    immediate_display = false;
                elseif p(2) <= 0 && p(2) <= p(1)
                    % insert in x
                    idx = find(p(1) >= x_thr, 1, 'last'); % never empty thanks to the -Inf
                    % x/y swap rather than insert?
                    if do_swap && strcmp(d_layout, 'y')
                        is_swap = ~mod(idx, 2);
                        idx = ceil(idx / 2);
                    else
                        is_swap = false;
                    end
                    if is_swap
                        new_layout_id.x = [layout_id_d.x(1:idx-1), dim_id, layout_id_d.x(idx+1:end)];
                        new_layout_id.y = [layout_id_d.y(1:didx-1), layout_id_d.x(idx), layout_id_d.y(didx:end)];
                    else
                        new_layout_id.x = [layout_id_d.x(1:idx-1), dim_id, layout_id_d.x(idx:end)];
                    end
                elseif p(1) <= 0 && p(2) <= .25
                    % goes in merged_data
                    if color_dim_id
                        % move away dimension that was occupying merged_data
                        % location
                        tmp = new_layout_id.(d_layout);
                        new_layout_id.(d_layout) = [tmp(1:didx-1), color_dim_id, tmp(didx:end)];
                        % and replace it by current dimension
                        new_layout_id.merged_data = dim_id;
                    else
                        % add current dimension to merged_data location
                        new_layout_id.merged_data(end+1) = dim_id;
                    end
                elseif p(1)<=0
                    % insert in y
                    idx = find(1 - p(2) >= y_thr, 1, 'last'); % never empty thanks to the +Inf
                    % x/y swap rather than insert?
                    if do_swap && strcmp(d_layout,'x')
                        is_swap = ~mod(idx, 2);
                        idx = ceil(idx / 2);
                    else
                        is_swap = false;
                    end
                    if is_swap
                        new_layout_id.x = [new_layout_id.x(1:didx-1), new_layout_id.y(idx), new_layout_id.x(didx:end)];
                        new_layout_id.y = [new_layout_id.y(1:idx-1), dim_id, new_layout_id.y(idx+1:end)];
                    else
                        new_layout_id.y = [new_layout_id.y(1:idx-1), dim_id, new_layout_id.y(idx:end)];
                    end
                elseif any(p>=1)
                    % xy and yx
                    if p(2) >= p(1)
                        % (zone above the graph)
                        [new_layout_id.xy, new_layout_id.yx] = deal(dim_id, []);
                    else
                        % (zone to the right of the graph)
                        [new_layout_id.xy, new_layout_id.yx] = deal([], dim_id);
                    end
                    if xy_dim_id
                        % move away dimension that was occupying xy or yx
                        % (swap with dim_id previous location)
                        tmp = new_layout_id.(d_layout);
                        new_layout_id.(d_layout) = [tmp(1:didx-1), xy_dim_id, tmp(didx:end)];
                    end
                else
                    % no change in layout, return
                    % note that there is no more zone for yx!, meaning that
                    % this display option is no longer accessible, it
                    % seemed to be too useless
                    return
                end
                % update organization (-> will trigger automatic display
                % update)
                if ~isequal(new_layout_id, layout_id)
                    layout_id = new_layout_id;
                    any_change = true;
                    L.D.set_layout_id(layout_id, immediate_display)
                end
                drawnow update
            end
            
            % finish 
            L.moving_dim = [];
            delete(L.moving_clone)
            L.moving_clone = [];
            if ~moved
                % label was not moved, make d the active x or y dimension
                % if appropriate
                make_dim_active(L.D, dim_id, 'toggle')
            elseif do_filter
                % apply the filter! -> this will cause a reslice and  a
                % global change in dimensions
                view_control.activate_inoperant_filter(dim_id)
            elseif isequal(layout_id, prev_layout_id)
                % do not put label back in its original place if it was only slightly moved: this allow
                % user to slightly move labels that cover other information
                if any_change
                    set_positions(L)
                end
            elseif L.do_immediate_display || isequal(layout_id,prev_layout_id)
                % update label positions once more to put the one for dimension
                % d in place
                set_positions(L)
            else
                % change organization (which will trigger data display
                % update)only now
                L.D.set_layout_id(layout_id,true); % second argument (true) is needed to force display update even though D.layout is already set to curlayout
            end
        end
    end
    
end
