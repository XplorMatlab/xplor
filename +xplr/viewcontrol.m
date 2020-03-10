classdef viewcontrol < xplr.graphnode
% view control
    
    properties (SetAccess='private')
        V               % parent 'view' object
        hp              % display panel
        items           % dimcontrols % uicontrols
        dimlist         % list of dimensions
        privatelists    % listcombo object
    end
    
    % Constructor
    methods
        function C = viewcontrol(V)
            % constructor viewcontrol
            
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
            
            % create a filter of key 1 for each dimension
            for i = 1:C.V.slicer.nddata
                C.dimaction('filter',1,i);
            end
        end
    end
    
    % Organization of items
    % items are organized vertically and are uicontrols or uipanels
    methods (Access='private')
        function init_items(C)
            fn_pixelsizelistener(C.hp,@(u,e)itemPositions(C))
            
            % note that other fields will be added, e.g. in addFilterItem
            C.items = struct('id',cell(1,0),'span',[],'obj',[]);
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
        function dimaction(C,flag,key,dim)
            % get list combo for the specified key
            isprivate = (key==0);
            
            % loop on selected dimensions
            newfilters = struct('d',cell(1,0),'F',[]);
            rmfilters = [];
            for d = dim
                % any filter already?
                slicefilters = C.V.slicer.filters;
                filteridx = find([slicefilters.dim]==d);
                if isempty(filteridx), filteridx = 0; end
                
%                 % list already shown?
%                 listidx = 0;
%                 if ~isempty(combo)
%                     shownfilters = combo.filters;
%                     for i=1:length(shownfilters)
%                         if isequal(getID(shownfilters(i).headerin),ID), listidx=i; break, end
%                     end
%                 end
                
                % create/replace/remove filter
                switch flag
                    case 'filter'
                        % create/replace filter

                        if filteridx
                            % a filter is already present in slicer, does
                            % it have the requested key?
                            F = slicefilters(filteridx).obj;
                            if F.linkkey == key
                                % if the existing filter has already the
                                % key of the requested new filter then
                                % nothing to do
                                continue
                            else
                                % the filter has to be remove before
                                % creating a new one with the new key
                                filterToRemove = C.V.slicer.filters(filteridx);
                                C.removefilter(filterToRemove);
                                % add the dimension to the slice's filters removal list
                                rmfilters = [rmfilters filteridx]; %#ok<AGROW>
                            end
                        end
                        % create the new filter
                        filterCreated = C.createfilter(d,key);
                        % add to the list of new filters
                        newfilters(end+1) = struct('d',d,'F',filterCreated); %#ok<AGROW>
                    case 'rmfilter'
                        if filteridx == 0
                            % no filter found for this dimension
                            continue
                        end
                        % remove filter from the viewcontrol and the bank
                        filter = C.V.slicer.filters(filteridx);
                        C.removefilter(filter);
                        
                         % add the dimension to the slice's filters removal list
                        rmfilters = [rmfilters filteridx]; %#ok<AGROW>
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
                    case 'showFilterWindow'
                        % re-create the filter window if not existing or
                        % put it in the front
                        
                        % if the filter is not private
                        if key ~= 0
                            % get filter from the filterSet
                            header = C.V.data.header(dim);
                            doshow=true;
                            filter = xplr.bank.getFilter(key,header,doshow);
                            xplr.bank.showList(filter);
                        end
                    case 'showFilterPointWindow'
                        % display Filter Point associated with this
                        % dimension
                        
                        if key ~= 0
                            % get filter from the filterSet
                            header = C.V.data.header(dim);
                            doshow=false;
                            F = xplr.bank.getFilter(key,header,doshow);
                            
                            % look for the filter with same key and with header
                            % corresponding wiht the headerout of the
                            % active filter
                            header = F.headerout;
                            doshow=true;
                            F = xplr.bank.getFilter(key,header,doshow);
                            % display the filter
                            xplr.bank.showList(F);
                        end
                end
            end
            
            % Add/remove all filters in the slicer at once (to have a unique
            % display update)
            if ~isempty(rmfilters)
                C.V.slicer.rmFilter(rmfilters)
            end
            if ~isempty(newfilters)
                C.V.slicer.addFilter({newfilters.d},[newfilters.F])
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
            
            backgroundColor = xplr.colors('linkkey',F.linkkey);
            panel.BackgroundColor = backgroundColor;
            % label
            hlab = uicontrol('parent',panel, ...
                'pos',[20 5 300 15], ...
                'style','text','string',label,'horizontalalignment','left', ...
                'backgroundcolor',backgroundColor, ...
                'enable', 'inactive', ...
                'buttondownfcn',@(u,e)clickFilterItem(C,d,id));
            C.items(itemidx).label = hlab;
           
            
            % buttons
            [ii jj] = ndgrid(-2:2);
            x = min(1,abs(abs(ii)-abs(jj))*.5);
            x(x==1) = NaN; x = repmat(x,[1 1 3]);
            
            % cross button to remove the filter
            rmFilterButton = uicontrol('parent',panel,'cdata',x, ...
                'unit', 'normalized', ...
                'position', [ 0.95 0.5 0.05 0.5 ], ...
                'callback',@(u,e)C.dimaction('rmfilter',1,d));
            fn_controlpositions(rmFilterButton, panel, [1 .5 0 .5], [-11 0 11 0]);
            
            
            % checkbox to disable and enable the filter
            uicontrol('parent',panel, ...
                'backgroundcolor',backgroundColor, ...
                'Style','checkbox', 'Value',1, ...
                'position', [ 6 6 13 12 ], ...
                'callback',@(u,e)C.dimaction('toggleactive',1,d));
            
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
            newidx = [];
            moved = fn_buttonmotion(@move,hf,'moved?');
            function move
                p = get(hf,'currentpoint'); p = p(1,2);
                newidx = fn_coerce( idx0 - round((p-p0)/ystep), 1, nfilter)
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
            
            % show filter if there was no move
            if ~moved, dimaction(C,'showFilterWindow',C.items(idxitem).F.linkkey,d), end
        end
    end
    
    % Private lists display
    methods (Access='private')
        function combo = getPrivateLists(C)
            combo = C.privatelists;
            controlorg = C.V.panels.allcontrols;
            % Create?
            if isempty(combo)
                disp 'warning: usage of private lists display has not been tested yet'
                combo = xplr.listcombo(C.V.panels.listcombo,0);
                C.privatelists = combo;
                connectlistener(combo,controlorg,'Empty',@(u,e)set(controlorg,'extents',[1 0]));
            end
            % Need to show it?
            if controlorg.extents(2) == 0
                % make combo visible
                controlorg.extents = [2 1];
            end
        end
    end
    
    % private filter management
    methods (Access='private')
        function removefilter(C,filter)
            % remove the filter from the viewcontrol and the bank
            % this function does not remove the filter from the slicer
            
            % if filter is empty, does nothing and leave the function
            if isempty(filter), return, end
            % remove filter from the items
            C.rmItem(['filter ' num2str(filter.dim)])
            % remove filter from the lists display
            % if the filter is private
            if filter.obj.linkkey == 0
                % remove the filter from the combo
                combo = C.getPrivateLists();
                combo.removeList(filter.obj)
            else
                % viewcontrol object C will be unregistered for the users
                % list of filter F; if this list will become empty, F will
                % be unregistered from the filters set
                xplr.bank.unregisterFilter(filter.obj,C) 
            end
        end
        
        function filter = createfilter(C,dimension,key)
            % create filter or get existing one from the
            % related public filters set
            
            header = C.V.data.header(dimension);
            % if the filter has to be private
            if key == 0
                % create private filter
                filter = xplr.filterAndPoint(header,'indices');
                % show filter in combo
                combo = C.getPrivateLists();
                combo.showList(filter)
            else
                % search for the filter in the bank with key and dimension
                doshow=true;
                filter = xplr.bank.getFilter(key,header,doshow,C);
            end
            
            % add the filter to the items, it is important that
            % filter.linkkey is set before using addFilterItem
            % TODO: change how the string filter is shifted
            str = ['filter ' header.label ' (' char(filter.F.slicefun) ')'];
            C.addFilterItem(dimension,str,filter)
        end
    end
    
    
end