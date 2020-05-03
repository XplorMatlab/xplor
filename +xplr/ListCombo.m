classdef ListCombo < hgsetget
    % function C = listcombo(container,filters)
    % if container is empty, a new figure will be created; this
    % figure will auto-delete once all filters will be removed
    
    properties
        container
        lists = xplr.List.empty(1, 0);
        filters = xplr.FilterAndPoint.empty(1, 0);
    end
    
    events
        Empty
    end
    
    % Constructor, add and remove lists
    methods
        function C = ListCombo(container, filters)
            
            % input
            if nargin < 1, container = []; end
            if nargin < 2, filters = []; end
                
            % create new containing figure? (in this case, set auto-delete)
            if isempty(container)
                screen_size = get(0, 'screen_size');
                container = figure('integerhandle', 'off', 'handlevisibility', 'off', ...
                    'numbertitle', 'off', 'menubar', 'none', ...
                    'name', 'Shared Filters', ...
                    'position', [min(80, screen_size(3)/20), max(screen_size(4)*.6-275, 45), 150, 550]);
                delete(findall(container, 'parent', container))
                add_listener(C, 'Empty', @(u,e)delete(container));
            end
            if ~isa(container, 'panelorganizer')
                container = panelorganizer(container, 'H');
                container.bordermode = 'push';
            end
            add_listener(container, 'ObjectBeingDestroyed', @(u,e)delete(C));
            C.container = container;
            
            % display lists
            C.add_list(filters)
        end
        function delete(C)
            if ~isprop(C, 'filters'), return, end
            if ~isempty(C.filters), notify(C, 'Empty'), end
        end
        function add_list(C, filter)
            % empty or multiple filters?
            if ~isscalar(filter)
                for i=1:length(filter), C.add_list(filter(i)), end
                return
            end
            
            % new panel
            [hp, id_x] = C.container.addSubPanel;
            
            % create list
            C.lists(id_x) = xplr.List(filter, 'in', hp);
            
            % remove panel when list will be distroyed
            add_listener(C.lists(id_x), 'ObjectBeingDestroyed', @(u,e)C.remove_list(hp));
            
            % memorize which filter is at this position
            C.filters(id_x) = filter;
        end
        function show_list(C, filter)
            % add filter list, only if not already present
            if fn_find(filter, C.filters)
                % filter is already shown
                return
            end
            add_list(C, filter)
        end
        function remove_list(C, x)
            % function remove_list(C,hp|id_x|filter)
            
            if ~isvalid(C), return, end
            
            if isa(x, 'matlab.ui.container.Panel') || isa(x, 'uipanel') || isnumeric(x)
                hp_or_id_x = x;
            elseif isa(x, 'xplr.DataOperand')
                hp_or_id_x = fn_find(x, C.filters);
                if isempty(hp_or_id_x)
                    return
                end
            else
                error argument
            end
            
            % remove panel
            id_x = C.container.remove_sub_panel(hp_or_id_x);
            C.lists(id_x) = [];
            C.filters(id_x) = [];
            % signal whether there are no more list being displayed
            if C.container.n_children == 0
                notify(C, 'Empty')
            end
        end
    end
    
    methods (Static)
        function C = test
            head = xplr.Header({'x', 10}, {'y', 12}, {'cond', {'a', 'b', 'c', 'd'}});
            for i=1:length(head)
                Filters(i) = xplr.FilterAndPoint(head(i)); %#ok<AGROW>
            end
            key = 2;
            C = xplr.ListCombo([], Filters);
            % repeat the 2d list
            C.add_list(Filters(2))
        end
    end
    
end
