classdef listcombo < hgsetget
    % function C = listcombo(container,filters)
    % if container is empty, a new figure will be created; this
    % figure will auto-delete once all filters will be removed
    
    properties
        container
        lists = xplr.list.empty(1,0);
        filters = xplr.filterAndPoint.empty(1,0);
    end
    
    events
        Empty
    end
    
    % Constructor, add and remove lists
    methods
        function C = listcombo(container,filters)
            
            % input
            if nargin<1, container = []; end
            if nargin<2, filters = []; end
                
            % create new containing figure? (in this case, set auto-delete)
            if isempty(container)
                screensize = get(0,'ScreenSize');
                container = figure('integerhandle','off','handlevisibility','off', ...
                    'numbertitle','off','menubar','none', ...
                    'name','Shared Filters', ...
                    'position',[min(80,screensize(3)/20) max(screensize(4)*.6-275,45) 150 550]);
                delete(findall(container,'parent',container))
                addlistener(C,'Empty',@(u,e)delete(container));
            end
            if ~isa(container,'panelorganizer')
                container = panelorganizer(container,'H');
                container.bordermode = 'push';
            end
            addlistener(container,'ObjectBeingDestroyed',@(u,e)delete(C));
            C.container = container;
            
            % display lists
            C.addList(filters)
        end
        function delete(C)
            if ~isprop(C,'filters'), return, end
            if ~isempty(C.filters), notify(C,'Empty'), end
        end
        function addList(C,filter)
            % empty or multiple filters?
            if ~isscalar(filter)
                for i=1:length(filter), C.addList(filter(i)), end
                return
            end
            
            % new panel
            [hp, idx] = C.container.addSubPanel;
            
            % create list
            C.lists(idx) = xplr.list(filter,'in',hp);
            
            % remove panel when list will be distroyed
            addlistener(C.lists(idx),'ObjectBeingDestroyed',@(u,e)C.removeList(hp));
            
            % memorize which filter is at this position
            C.filters(idx) = filter;
        end
        function showList(C,filter)
            % add filter list, only if not already present
            if fn_find(filter, C.filters)
                % filter is already shown
                return
            end
            addList(C,filter)
        end
        function removeList(C,x)
            % function removeList(C,hp|idx|filter)
            
            if ~isvalid(C), return, end
            
            if isa(x,'matlab.ui.container.Panel') || isa(x,'uipanel') || isnumeric(x)
                hp_or_idx = x;
            elseif isa(x,'xplr.dataOperand')
                hp_or_idx = fn_find(x,C.filters);
                if isempty(hp_or_idx)
                    return
                end
            else
                error argument
            end
            
            % remove panel
            idx = C.container.removeSubPanel(hp_or_idx);
            C.lists(idx) = [];
            C.filters(idx) = [];
            % signal whether there are no more list being displayed
            if C.container.nchildren==0
                notify(C,'Empty')
            end
        end
    end
    
    methods (Static)
        function C = test
            head = xplr.header({'x' 10},{'y' 12},{'cond' {'a' 'b' 'c' 'd'}});
            for i=1:length(head)
                Filters(i) = xplr.filterAndPoint(head(i)); %#ok<AGROW>
            end
            key = 2;
            C = xplr.listcombo([],Filters);
            % repeat the 2d list
            C.addList(Filters(2))
        end
    end
    
end


