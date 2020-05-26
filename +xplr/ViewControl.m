classdef ViewControl < xplr.GraphNode
% view control
    
    properties (SetAccess='private')
        V               % parent 'view' object
        hp              % display panel
        items           % dimcontrols % uicontrols
        dim_list         % list of dimensions
        private_lists    % list_combo object
        context_menu         % context menu
    end
    
    % Constructor
    methods
        function C = ViewControl(V)
            % constructor viewcontrol
            
            % parent 'view' object and panel
            C.V = V;
            C.hp = V.panels.control;
            
            % items
            init_items(C)
            % (data)
            data_str = 'Data';
            if ~isempty(V.data.name), data_str = [data_str, ' (', V.data.name, ')']; end
            C.new_item('data', 1, ...
                {'style', 'text', 'string', V.data.name, ...
                'backgroundcolor', xplr.colors('gui.controls.dataname'), ...
                'enable', 'inactive', 'buttondownfcn', @(u,e)C.data_context_menu()})
            % (list of data dimensions)
            C.dim_list = C.new_item('dim_list', 4, ...
                {'style', 'listbox', 'string', {V.data.header.label}, 'max', 2, ...
                'callback', @(u,e)C.dimension_context_menu()});
            C.context_menu = uicontextmenu(V.hf);
            
            % some changes needed when data header is changed
            C.add_listener(V.data, 'ChangedData', @(u,e)data_change(C, e));
            
            % create initial list of filters
            % (determine which filters should be active for the slice to be
            % displayable)
            nd = C.V.data.nd;
            active = false(1, nd);
            n_dim_max = min(4, nd); % no more than 4 dimensions visible
            active(n_dim_max+1:end) = true;
            for i = n_dim_max:-1:1
                % test displayable
                sz = C.V.data.sz;
                sz(active) = 1;
                display_mode = C.V.D.display_mode;
                layout = xplr.DisplayLayout(C.V.D).dimension_number();
                if xplr.ViewDisplay.test_displayable(sz, display_mode, layout)
                    break
                end
                % not displayable -> activate one filter more, starting
                % from the end
                active(i) = true;
            end
            % (add filters)
            key = 1;
            if any(active)
                C.dim_action('addfilter', num2cell(find(active)), key) 
            end
        end
    end
    
    % Organization of items
    % items are organized vertically and are uicontrols or uipanels
    methods (Access='private')
        function init_items(C)
            fn_pixelsizelistener(C.hp, @(u,e)item_positions(C))
            
            % note that other fields will be added, e.g. in add_filter_item
            C.items = struct('id', cell(1, 0), 'span', [], 'obj', []);
        end
        function item_positions(C, idx)
            if nargin<2, idx = 1:length(C.items); end
            [W, H] = fn_pixelsize(C.hp);
            h = 22; % item height, in pixel
            dx = 2;
            dy = 2;
            w_max = Inf;
            w = max(1, min(w_max, W-2*dx));
            x0 = (W-w)/2;
            y_starts = [0, cumsum([C.items.span])];
            for i=row(idx)
                yspan = C.items(i).span;
                set(C.items(i).obj, 'units', 'pixel', 'position', [x0, H-(y_starts(i)+yspan)*(h+dy), w, yspan*h+(yspan-1)*dy])
            end
        end
        function [obj, idx] = new_item(C, id, span, control_prop)
            % function [obj idx] = new_item(C,id,span[,{uicontrol properties}])
            % function [obj idx] = new_item(C,id,span,'panel')
            if nargin<4 || iscell(control_prop)
                if nargin<4, control_prop = {}; end
                obj = uicontrol('parent', C.hp, ...
                    'backgroundcolor', xplr.colors('gui.controls.item'), ...
                    control_prop{:});
            elseif strcmp(control_prop, 'panel')
                obj = uipanel('parent', C.hp, 'bordertype', 'none', 'units', 'pixel');
            end
            idx = length(C.items) + 1;
            C.items(idx).pos = sum([C.items.span]) + 1;
            C.items(idx).id = id;
            C.items(idx).span = span;
            C.items(idx).obj = obj;
            item_positions(C, idx)
            if nargout == 0, clear obj, end
        end
        function rm_item(C, id)
            idx = fn_find(id, {C.items.id});
            delete_valid([C.items(idx).obj])
            C.items(idx) = [];
            item_positions(C)
        end
    end
    
    % Data (edit headers)
    methods
        function data_context_menu(C)
            % init context menu
            m = C.context_menu;
            delete(get(m, 'children'))
            
            % create entries
            uimenu(m, 'label', 'Edit header information', ...
                'callback', @(u,e)C.edit_header())
            uimenu(m, 'label', 'Open data in a new xplor window', ...
                'callback', @(u,e)xplor(C.V.data))
            
            % make menu visible
            p = get(C.V.hf, 'currentpoint');
            p = p(1, 1:2);
            set(m, 'Position', p, 'Visible', 'on')
        end
        function edit_header(C)
            data = C.V.data;
            cur_head = data.header;
            new_head = xplr.EditHeader(C.V.data);
            if isempty(new_head), return, end % user closed window: cancel
            dim_chg = false(1,data.nd);
            for i=1:data.nd, dim_chg(i) = ~isequal(new_head(i), cur_head(i)); end
            if any(dim_chg)
                dim = find(dim_chg);
                new_head_dim = xplr.DimHeader(new_head(dim));
                C.V.data.update_data('chg_dim', dim, [], data.data, new_head_dim)
            end
        end
        function data_change(C, e)
            switch e.flag
                case 'global'
                    error 'global data change not handled'
                case 'chg_dim'
                    % update dimension list
                    set(C.dim_list, 'string', {C.V.data.header.label})
                otherwise
                    % no change needed
            end
        end
    end
    
    % Dimensions menu
    methods
        function dimension_context_menu(C)
            % init context menu
            m = C.context_menu;
            delete(get(m, 'children'))
            
            % selected dimension(s)
            dim = get(C.dim_list, 'value');
            dim_id = [C.V.data.header(dim).dim_id];
            
            % some strings to handle singular vs. plural
            dim_str = fn_switch(isscalar(dim_id), 'this dimension', 'these dimensions');
            filter_str = fn_switch(isscalar(dim_id), 'filter', 'filters');
            
            % add or change filter(s)
            % (2D with key 1)
            if length(dim_id) == 2
                % (using key 1)
                uimenu(m, 'label', 'Filter with shared 2D filter', 'separator', 'on', ...
                    'callback', @(u,e)dim_action(C, 'addfilter', {dim_id}, 1))
                next_separator = 'off';
            else
                next_separator = 'off';
            end
            % (1D with key 1)
            uimenu(m, ...
                'label', ['Filter with shared 1D ', filter_str], 'separator', next_separator, ...
                'callback', @(u,e)dim_action(C, 'addfilter', num2cell(dim_id), 1))
            % (more options: select among available keys)
            available_keys = xplr.Bank.available_filter_keys('filterAndPoint');
            new_key = max(available_keys) + 1;
            key_values = [setdiff(available_keys, 1), new_key];
            m2 = uimenu(m, 'label', 'Filter with');
            if length(dim_id) == 2
                for key_value = key_values
                    uimenu(m2, 'label', ['shared 2D filter ', num2str(key_value)], ...
                        'callback', @(u,e)dim_action(C, 'addfilter', {dim_id}, key_value));
                end
            end
            for key_value = key_values
                uimenu(m2, 'label', ['shared 1D ', filter_str, ' ', num2str(key_value)], ...
                    'callback', @(u,e)dim_action(C, 'addfilter', num2cell(dim_id), key_value));
            end
            uimenu(m2, 'label', ['private 1D ', filter_str], ...
                'callback', @(u,e)dim_action(C, 'addfilter', num2cell(dim_id), 0))
            
            % remove filters in these dimensions
            uimenu(m, 'label', ['Remove ', filter_str], 'separator', 'on', ...
                'callback', @(u,e)dim_action(C, 'rmfilter', dim_id))
            
            % filter all others dimension
            uimenu(m, ...
                'label', ['View ', dim_str, ', filter others'], 'separator', 'on', ...
                'callback', @(u,e)dim_action(C, 'viewdim', dim_id, 1))
            uimenu(m, 'label', ['View ', dim_str, ' in a new window'], ...
                'callback', @(u,e)dim_action(C, 'newwindow_viewdim', dim_id, 1))
            
            % make menu visible
            p = get(C.V.hf, 'currentpoint');
            p = p(1, 1:2);
            set(m, 'Position', p, 'Visible', 'on')
        end
        function dim_action(C, flag, dim_id, varargin)
            % function dim_action(C,'addfilter',dim_ids[,key[,active]])
            % function dim_action(C,'rmfilter|showfilter',dim_id)
            % function dim_action(C,'set_active',dim_id,value)
            %---
            % if flag is 'addfilter', dims can be a cell array, to defined
            % several filters at once for example
            % dim_action(C,'addfilter',{[1 2] 3}) will add two filters,
            % first a 2D filter in dimensions [1 2], second a 1D filter in
            % dimension 3
            %
            % dim_id is supposed to be the unique identifier of some
            % dimension(s), but for commodity it can also be the dimension
            % number, or the dimension label
            
            % other window
            if strfind(flag, 'newwindow') %#ok<STRIFCND>
                % open data in a new window: flag can be either
                % 'otherwindow' or 'otherwindow_action' where 'action' is
                % to be executed in this window
                V2 = xplor(C.V.data);
                tokens = regexp(flag, 'newwindow_(.*)', 'tokens');
                if ~isempty(tokens)
                    V2.C.dim_action(tokens{1}{1}, dim_id, varargin{:})
                end
                return
            end
            
            % convert dimension numbers or labels to dimension identifiers
            dim_id = C.V.data.dimension_id(dim_id);
            
            % 'addfilter' flag -> several filters at once
            if strcmp(flag, 'addfilter')
                % dims will be a cell array: list of dimensions, per filter
                % dim_id will be an array: list of all affected dimensions
                if ~iscell(dim_id)
                    if ~isscalar(dim_id), error 'array of dim_id values is ambiguous, use a cell array instead', end
                    dim_ids = {dim_id};
                else
                    dim_ids = dim_id; % several set of one or several dimensions
                    dim_id = unique([dim_ids{:}]);
                end
                if length(dim_id) < length([dim_ids{:}])
                    error 'some dimension is repeated in filter(s) definition'
                end
            end
            
            % list of filters in the selected dimensions
            filters_idx = find(fn_map({C.V.slicer.filters.dim_id}, @(dd)any(ismember(dd, dim_id)), 'array'));
            current_filters_dim = C.V.slicer.filters(filters_idx); % current filters acting on dimensions within dd

            % filters to remove
            if ismember(flag, {'addfilter', 'rmfilter', 'viewdim'})
                % remove filter from the viewcontrol and the bank
                for filter = current_filters_dim
                    C.remove_filter(filter);
                end
                
                % remove filters from the slicer
                do_slicing = strcmp(flag, 'rmfilter'); % no need to reslice yet for 'addfilter', reslice will occur when adding the new filter(s)
                C.V.slicer.rm_filter(filters_idx, do_slicing);
            end
            
            % filters to add
            if ismember(flag, {'addfilter', 'viewdim'})
                if nargin >= 4, key = varargin{1}; else, key = 1; end
                if nargin >= 5, active = varargin{2}; else, active = true; end
                if strcmp(flag, 'addfilter')
                    dim_ids_add = dim_ids; % already a cell array
                else
                    % add 1D filters il all dimensions that we do not want
                    % to view and that are not already filtered
                    no_view_dim_id = setdiff([C.V.data.header.dim_id], dim_id, 'stable');
                    cur_filt_dim_id = [C.V.slicer.filters.dim_id];
                    dim_ids_add = setdiff(no_view_dim_id, cur_filt_dim_id, 'stable');
                    % among these dimensions, attempt to find pairs of
                    % measure headers with same units to set 2D filter
                    % instead of two 1D filters
                    head = C.V.data.header_by_id(dim_ids_add);
                    connections = measure_grouping(head);
                    pairs = {};
                    while any(connections(:))
                        [i, j] = find(connections, 1, 'first');
                        pairs{end+1} = dim_ids_add(sort([i, j]));
                        connections([i, j], :) = false;
                        connections(:, [i, j]) = false;
                    end
                    dim_ids_add = [pairs, num2cell(setdiff(dim_ids_add, [pairs{:}], 'stable'))];
                end
                n_add = length(dim_ids_add);
                if n_add > 0
                    if n_add>1 && isscalar(key), key = repmat(key, 1, n_add); end
                    if n_add>1 && isscalar(active), active = repmat(active, 1, n_add); end
                    % loop on dimension sets
                    new_filters = struct('dim_id', cell(1, 0), 'F', [], 'active', []);
                    for i = 1:length(dim_ids_add)
                        F = C.create_filter_and_item(dim_ids_add{i}, key(i), active(i));
                        new_filters(end+1) = struct('dim_id', dim_ids_add{i}, 'F', F, 'active', active(i)); %#ok<AGROW>
                    end
                    C.V.slicer.add_filter({new_filters.dim_id}, [new_filters.F], [new_filters.active]) % slicing will occur now
                else
                    % we might have removed filters before without updating
                    % completely the slice
                    C.V.slicer.apply_pending()
                end
                
                % adjust display mode and layout if it seems appropriate
                D = C.V.D;
                if strcmp(flag, 'viewdim')
                    if isscalar(dim_id)
                        D.set_dim_location(dim_id, 'x', strcmp(D.display_mode, 'time courses'))
                        D.display_mode = 'time courses';
                    elseif length(dim_id) == 2
                        D.set_dim_location(dim_id, {'x', 'y'}, strcmp(D.display_mode, 'image'))
                        D.display_mode = 'image';
                    end
                else
                    nsdim_id = non_singleton_dim_id(C.V.slice.header);
                    if isscalar(nsdim_id)
                        D.set_dim_location(nsdim_id, 'x', strcmp(D.display_mode, 'time courses'))
                        D.display_mode = 'time courses';
                    end
                end
            end
            
            % show filter, set filter active
            switch flag
                case 'set_active'
                    active = varargin{1};
                    % show label(s) as enabled/disabled
                    for filter = current_filters_dim
                        item_idx = fn_find({'filter', filter.dim_id}, {C.items.id});
                        h_lab = C.items(item_idx).label;
                        set(h_lab, 'enable', fn_switch(active, 'inactive', 'off'))
                        drawnow
                    end
                    % toggle filter active in slicer
                    C.V.slicer.chg_filter_active(filters_idx, active)
                case 'showfilter'
                    for filter = current_filters_dim
                        F = filter.obj;
                        if ~isscalar(filter.dim_id)
                            disp('cannot display list for ND filter')
                        elseif F.link_key == 0
                            % private filter
                            combo = C.get_private_lists();
                            combo.show_list(F)
                        else
                            xplr.Bank.show_list(F);
                        end
                    end
            end

            % Empty the dimension selection
            set(C.dim_list, 'value', [])
        end
        function add_filter_item(C, dim_id, F, active)
            % panel
            id = {'filter', dim_id};
            [panel, item_idx] = C.new_item(id, 1, 'panel');
            background_color = xplr.colors('link_key', F.link_key);
            panel.BackgroundColor = background_color;
            
            % store the filter
            C.items(item_idx).F = F;
            
            % label
            label = ['filter ', F.header_in.label, ' (', F.F.slice_fun_str, ')'];
            h_lab = uicontrol('parent', panel, ...
                'position', [20, 5, 300, 15], ...
                'style', 'text', 'string', label, 'horizontalalignment', 'left', ...
                'backgroundcolor', background_color, ...
                'enable', fn_switch(active, 'inactive', 'off'), ...
                'buttondownfcn', @(u,e)click_filter_item(C, dim_id, id), ...
                'uicontextmenu', uicontextmenu(C.V.hf, 'callback', @(m, e)F.context_menu(m)));
            C.items(item_idx).label = h_lab;
            
            % change label upon operation change
            function check_operation_change(~, e)
                if strcmp(e.type, 'operation')
                    label = ['filter ', F.header_in.label, ' (', F.F.slice_fun_str, ')'];
                    set(h_lab, 'string', label)
                end
            end
            connect_listener(F.F,h_lab, 'ChangedOperation', @check_operation_change);
            
            % buttons
            [ii, jj] = ndgrid(-2:2);
            x = min(1, abs(abs(ii) - abs(jj))*.5);
            x(x == 1) = NaN;
            x = repmat(x, [1, 1, 3]);
            
            % cross button to remove the filter
            rm_filter_button = uicontrol('parent', panel, 'cdata', x, ...
                'unit', 'normalized', ...
                'position', [0.95, 0.5, 0.05, 0.5], ...
                'callback', @(u,e)C.dim_action('rmfilter', dim_id));
            fn_controlpositions(rm_filter_button, panel, [1, .5, 0, .5], [-11, 0, 11, 0]);
            
            % checkbox to disable and enable the filter
            uicontrol('parent', panel, ...
                'backgroundcolor', background_color, ...
                'Style', 'checkbox', 'Value', active, ...
                'position', [6, 6, 13, 12], ...
                'callback', @(u,e)C.dim_action('set_active', dim_id, get(u, 'value')));
            
        end
        function click_filter_item(C, d, id)
            hf = C.V.hf;
            switch get(hf, 'selectiontype')
                case 'normal'
                    % try to move the filter, if no move, toggle active:
                    % see the code later
                otherwise
                    return
            end
                                
            % get items corresponding to filters
            idx_filter = find(~fn_isemptyc({C.items.F}));
            if ~all(diff(idx_filter) == 1), error 'filters should be contiguous', end
            filter_items = C.items(idx_filter);
            n_filter = length(idx_filter);
            
            % index and position of selected filter
            idx_item = fn_find(id, {C.items.id});
            idx_0 = idx_item - (idx_filter(1) - 1);
            idx_other = setdiff(1:n_filter, idx_0);
            obj = C.items(idx_item).obj;
            pos_0 = get(obj, 'position');
            y_step = 24;
            
            % move
            p0 = get(hf, 'currentpoint');
            p0 = p0(1, 2); % only vertical position matters
            new_idx = [];
            moved = fn_buttonmotion(@move, hf, 'moved?');
            function move
                p = get(hf, 'currentpoint');
                p = p(1, 2);
                new_idx = fn_coerce(idx_0 - round((p-p0)/y_step), 1, n_filter);
                % set all items position
                C.items(idx_filter) = filter_items([idx_other(1:new_idx-1), idx_0 idx_other(new_idx:end)]);
                C.item_positions
                % set selected item position
                new_pos = pos_0;
                new_pos(2) = pos_0(2) + fn_coerce(p-p0, [idx_0-n_filter idx_0-1]*y_step);
                set(obj, 'position', new_pos)
            end
            if moved
                % re-position correctly the selected item
                C.item_positions
                % apply filters permutation
                perm = [idx_other(1:new_idx-1), idx_0, idx_other(new_idx:end)];
                C.V.slicer.perm_filters(perm)
            end
            
            % show filter if there was no move
            if ~moved, dim_action(C, 'showfilter', d), end
        end
    end
    
    % Filters display
    methods (Access='private')
        function remove_filter(C, filter)
            % remove the filter from the viewcontrol and the bank
            % this function does not remove the filter from the slicer
            
            % if filter is empty, does nothing and leave the function
            if isempty(filter), return, end
            % remove filter from the items
            C.rm_item({'filter', filter.dim_id})
            % remove filter from the lists display
            % if the filter is private
            if filter.obj.link_key == 0
                % remove the filter from the combo
                combo = C.get_private_lists();
                combo.remove_list(filter.obj)
            else
                % viewcontrol object C will be unregistered for the users
                % list of filter F; if this list will become empty, F will
                % be unregistered from the filters set
                xplr.Bank.unregister_filter(filter.obj, C)
            end
        end
        function F = create_filter_and_item(C, dim_id, key, active)
            % create filter or get existing one from the
            % related public filters set
            
            header = C.V.data.header_by_id(dim_id);
            % if the filter has to be private
            if key == 0
                % create private filter
                F = xplr.FilterAndPoint(header);
                % show filter in combo
                if isscalar(dim_id)
                    combo = C.get_private_lists();
                    if active, combo.show_list(F), end
                end
            else
                % search for the filter in the bank with key and dimension
                do_show = true;
                F = xplr.Bank.get_filter_and_point(key, header, C, do_show);
            end
            
            % add the filter to the items, it is important that
            % filter.link_key is set before using add_filter_item
            % TODO: change how the string filter is shifted
            C.add_filter_item(dim_id, F, active)
        end
    end
    
    % Private lists display
    methods (Access='private')
        function combo = get_private_lists(C)
            combo = C.private_lists;
            control_org = C.V.panels.all_controls;
            % Create?
            if isempty(combo)
                disp 'warning: usage of private lists display has not been tested yet'
                combo = xplr.ListCombo(C.V.panels.list_combo);
                C.private_lists = combo;
                connect_listener(combo, control_org, 'Empty', @(u,e)set(control_org, 'extents', [1, 0]));
            end
            % Need to show it?
            if control_org.extents(2) == 0
                % make combo visible
                control_org.extents = [2, 1];
            end
        end
    end
    
    
end
