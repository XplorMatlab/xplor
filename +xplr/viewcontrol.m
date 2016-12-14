classdef viewcontrol < hgsetget
    
    properties (SetAccess='private')
        V   % parent 'view' object
        hp  % display panel
        items
        %         dimcontrols % uicontrols
        dimlist     % list of dimensions
        lists       % listcombo object
    end
    
    % Constructor
    methods
        function C = viewcontrol(V)
            % parent 'view' object and panel
            C.V = V;
            C.hp = V.panels.control;
            
            % items
            init_items(C)
            % (data)
            C.newItem('data',1,{'string',V.data.name,'backgroundcolor',xplr.colors('gui.controls.dataname'), ...
                'callback',@(u,e)editHeader(C)})
            % (list of data dimensions)
            C.dimlist = C.newItem('dimlist',4,{'style','listbox','string',{V.data.header.label},'max',2, ...
                'callback',@(u,e)C.V.context.raise('datadim',get(u,'value'))});            
        end
    end
    
    % Organization of items
    % items are organized vertically and are uicontrols or uipanels
    methods (Access='private')
        function init_items(C)
            fn_pixelsizelistener(C.hp,@(u,e)itemPositions(C))
            C.items = struct('id',cell(1,0),'span',[],'obj',[]); % note that other fields will be added, e.g. in addFilterItem
        end
        function itemPositions(C,idx)
            if nargin<2, idx = 1:length(C.items); end
            [W H] = fn_pixelsize(C.hp);
            h = 22; % item height, in pixel
            dx = 2; dy = 2;
            wmax = Inf;
            w = max(1,min(wmax,W-2*dx));
            x0 = (W-w)/2;
            ystarts = [0 cumsum([C.items.span])];
            for i=row(idx)
                yspan = C.items(i).span;
                set(C.items(i).obj,'units','pixel','pos',[x0 H-(ystarts(i)+yspan)*(h+dy) w yspan*h+(yspan-1)*dy])
            end
        end
        function [obj idx] = newItem(C,id,span,controlprop)
            % function [obj idx] = newItem(C,id,span[,{uicontrol properties}])
            % function [obj idx] = newItem(C,id,span,'panel')
            if nargin<4 || iscell(controlprop)
                obj = uicontrol('parent',C.hp, ...
                    'backgroundcolor',xplr.colors('gui.controls.item'), ...
                    controlprop{:});
            elseif strcmp(controlprop,'panel')
                obj = uipanel('parent',C.hp,'bordertype','none','units','pixel');
            end
            idx = length(C.items)+1;
            C.items(idx).pos = sum([C.items.span])+1;
            C.items(idx).id = id;
            C.items(idx).span = span;
            C.items(idx).obj = obj;
            itemPositions(C,idx)
            if nargout==0, clear obj, end
        end
        function rmItem(C,id)
            idx = strcmp({C.items.id},id);
            deleteValid([C.items(idx).obj])
            C.items(idx) = [];
            itemPositions(C)
        end
    end
    
    % Data (edit headers)
    methods
        function editHeader(C)
            data = C.V.data;
            curhead = data.header;
            newhead = xplr.editHeader(C.V.data);
            if isempty(newhead), return, end % user closed window: cancel
            dimchg = false(1,data.nd);
            for i=1:data.nd, dimchg(i) = ~isequal(newhead(i),curhead(i)); end
            if any(dimchg)
                dim = find(dimchg);
                C.V.data.updateData('chgdim',dim,[],data.data,newhead(dim))
            end
        end
    end
    
    % Dimensions menu
    methods
        function dimaction(C,flag,dim)
            % init lists if needed
            if isempty(C.lists) && fn_ismemberstr(flag,'filter')
                init_lists(C)
            end
            
            % loop on selected dimensions
            newfilters = struct('d',cell(1,0),'F',[]);
            rmfilters = [];
            for d = dim
                head = C.V.data.header(d);
                
                % any filter already?
                slicefilters = C.V.slicer.filters;
                filteridx = find([slicefilters.dim]==d);
                if isempty(filteridx), filteridx = 0; end
                
                % list already shown?
                listidx = 0;
                if ~isempty(C.lists)
                    shownfilters = C.lists.filters;
                    for i=1:length(shownfilters)
                        if shownfilters(i).shared.privatedim==d, listidx=i; break, end
                    end
                end
                
                % add/show/remove filter
                switch flag
                    case 'filter'
                        % create/get filter
                        if ~filteridx
                            if listidx
                                if fn_dodebug, warning 'list should not be present if the data is not filtered in this dimension', end
                                F = shownfilters(listidx);
                            else                            % create new filter and a list acting on it
                                F = xplr.filterAndPoint(head,'indices');
                                % add information about which dim is being filtered
                                F.shared.isprivate = true;
                                F.shared.privatedim = d;
                            end
                            % add to the list of new filters
                            newfilters(end+1) = struct('d',d,'F',F); %#ok<AGROW>
                            % add the filter to the items
                            str = ['filter ' head.label ' (' char(F.F.slicefun) ')'];
                            C.addFilterItem(d,str,F)
                        else
                            F = slicefilters(filteridx).obj;
                        end
                        % add the filter to the lists
                        if ~listidx
                            C.lists.addList(F)
                        end
                        figure(C.lists.container.hobj) % show the figure with filters
                    case 'rmfilter'
                        % add the dimension to the filter removal list
                        if filteridx, rmfilters = [rmfilters filteridx]; end %#ok<AGROW>
                        % remove filter from the items
                        C.rmItem(['filter ' num2str(d)])
                        % remove filter from the lists
                        if listidx, C.lists.removeList(listidx), end
                    case 'toggleactive'
                        itemidx = find(strcmp(['filter ' num2str(d)],{C.items.id}));
                        hlab = C.items(itemidx).label;
                        % active or inactive?
                        active = strcmp(get(hlab,'enable'),'off');
                        % enable/disable label
                        set(hlab,'enable',fn_switch(active,'inactive','off'))
                        % toggle filter active in slicer
                        disp 'very bad way to determine filter index'
                        filteridx = itemidx-2;
                        C.V.slicer.chgFilterActive(filteridx,active)
                end
            end
            
            % Add/remove all filters in the slicer at once (to have a unique
            % display update)
            if ~isempty(newfilters)
                C.V.slicer.addFilter({newfilters.d},[newfilters.F])
            end
            if ~isempty(rmfilters)
                C.V.slicer.rmFilter(rmfilters)
            end
            
            % Empty the dimension selection
            set(C.dimlist,'value',[])
        end
        function addFilterItem(C,d,label,F)
            % panel
            id = ['filter ' num2str(d)];
            [panel itemidx] = C.newItem(id,1,'panel');
            % store the filter
            C.items(itemidx).F = F;
            % label
            hlab = uicontrol('parent',panel,'units','normalized','pos',[0 0 1 1], ...
                'style','text','string',label,'horizontalalignment','left', ...
                'backgroundcolor',xplr.colors('linkkey',1), ...
                'enable','inactive','buttondownfcn',@(u,e)clickFilterItem(C,d,id));
            C.items(itemidx).label = hlab;
            % buttons
            [ii jj] = ndgrid(-2:2);
            x = min(1,abs(abs(ii)-abs(jj))*.5);
            x(x==1) = NaN; x = repmat(x,[1 1 3]);
            uclose = uicontrol('parent',panel,'cdata',x, ...
                'callback',@(u,e)C.dimaction('rmfilter',d));
            fn_controlpositions(uclose,panel,[1 .5 0 .5],[-11 0 11 0])
        end
        function clickFilterItem(C,d,id)
            hf = C.V.hf;
            switch get(hf,'selectiontype')
                case 'normal'
                    % try to move the filter, if no move, toggle active:
                    % see the code later
                otherwise
                    return
            end
                                
            % get items corresponding to filters
            idxfilter = find(~fn_isemptyc({C.items.F}));
            if ~all(diff(idxfilter)==1), error 'filters should be contiguous', end
            filteritems = C.items(idxfilter);
            nfilter = length(idxfilter);
            
            % index and position of selected filter
            idxitem = find(strcmp({C.items.id},id));
            idx0 = idxitem-(idxfilter(1)-1);
            idxother = setdiff(1:nfilter,idx0);
            obj = C.items(idxitem).obj;
            pos0 = get(obj,'pos');
            ystep = 24;
            
            % move
            p0 = get(hf,'currentpoint'); p0 = p0(1,2); % only vertical position matters
            moved = fn_buttonmotion(@move,hf,'moved?');
            function move
                p = get(hf,'currentpoint'); p = p(1,2);
                newidx = fn_coerce( idx0 - round((p-p0)/ystep), 1, nfilter);
                % set all items position
                C.items(idxfilter) = filteritems([idxother(1:newidx-1) idx0 idxother(newidx:end)]);
                C.itemPositions
                % set selected item position
                newpos = pos0; newpos(2) = pos0(2) + fn_coerce(p-p0,[idx0-nfilter idx0-1]*ystep);
                set(obj,'pos',newpos)
            end
            if moved
                % re-position correctly the selected item
                C.itemPositions
                % apply filters permutation
                perm = [idxother(1:newidx-1) idx0 idxother(newidx:end)];
                C.V.slicer.permFilters(perm)
            end
            
            % toggle active if there was no move
            if ~moved, dimaction(C,'toggleactive',d), end
        end
    end
    
    % Lists
    methods (Access='private')
        function init_lists(C)
            C.lists = xplr.listcombo();
            addlistener(C.lists,'Empty',@(u,e)set(C,'lists',[]));
            %             % local lists
            %             C.lists = xplr.listcombo(C.V.panels.listcombo,0);
            %             controlorg = C.V.panels.allcontrols;
            %             addlistener(C.lists,'Empty',@(u,e)set(controlorg,'extents',[1 0]));
        end
    end
    
end