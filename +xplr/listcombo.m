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
            [hp idx] = C.container.addSubPanel;
            
            % create graphic objects
            % (list)
            hlist = uicontrol('parent',hp,'style','listbox');
            fn_controlpositions(hlist,hp,[0 0 1 1],[8 5 -16 -5-21-2])
            % (label)
            hlabel = uicontrol('parent',hp,'style','text', ...
                'backgroundcolor',xplr.colors('linkkey',filter.linkkey));
            fn_controlpositions(hlabel,hp,[0 1 1 0],[8 -21 -8-18 18])
            % (close button)
            x = fn_printnumber(ones(18),'x','pos','center')'; 
            x(x==1) = NaN; x = repmat(x,[1 1 3]);
            hclose = uicontrol('parent',hp,'cdata',x,'callback',@(u,e)C.removeList(hp));
            fn_controlpositions(hclose,hp,[1 1],[-8-18 -3-18 18 18])

            % create list
            C.lists(idx) = xplr.list(filter,'in',[hlist hlabel]);
            
            % memorize which filter is at this position
            C.filters(idx) = filter;
            
            % watch filter deletion
            addlistener(filter,'ObjectBeingDestroyed',@(u,e)removeList(C,hp));
        end
        function showList(C,filter)
            % add filter list, only if not already present
            hID = getID(filter.headerin);
            for idx = 1:length(C.filters)
                if isequal(hID,getID(C.filters(idx).headerin)), return, end
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


