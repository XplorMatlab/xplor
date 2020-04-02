classdef viewcontrol < xplr.graphnode
% view control
    
    properties (SetAccess='private')
        V               % parent 'view' object
        hp              % display panel
        items           % dimcontrols % uicontrols
        dimlist         % list of dimensions
        privatelists    % listcombo object
        contextmenu         % context menu
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
            C.newItem('data',1, ...
                {'style','text','string',V.data.name, ...
                'backgroundcolor',xplr.colors('gui.controls.dataname'), ...
                'enable','inactive','buttondownfcn',@(u,e)C.dataContextMenu()})
            % (list of data dimensions)
            C.dimlist = C.newItem('dimlist',4, ...
                {'style','listbox','string',{V.data.header.label},'max',2, ...
                'callback',@(u,e)C.dimensionContextMenu()});
            C.contextmenu = uicontextmenu(V.hf);
            
            % some changes needed when data header is changed
            C.addListener(V.data,'ChangedData',@(u,e)datachange(C,e));            
            
            % create initial list of filters
            % (determine which filters should be active for the slice to be
            % displayable)
            nd = C.V.data.nd;
            active = false(1,nd);
            ndimmax = min(4,nd); % no more than 4 dimensions visible
            active(ndimmax+1:end) = true;
            for i = ndimmax:-1:1
                % test displayable
                sz = C.V.data.sz; sz(active) = 1;
                displaymode = C.V.D.displaymode;
                layout = xplr.displaylayout(C.V.D).dimensionNumber();
                if xplr.viewdisplay.testDisplayable(sz,displaymode,layout)
                    break
                end
                % not displayable -> activate one filter more, starting
                % from the end
                active(i) = true;
            end
            % (add filters)
            key = 1;
            if any(active)
                C.dimaction('addfilter',num2cell(find(active)),key) 
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
                set(C.items(i).obj,'units','pixel','position',[x0 H-(ystarts(i)+yspan)*(h+dy) w yspan*h+(yspan-1)*dy])
            end
        end
        function [obj, idx] = newItem(C,id,span,controlprop)
            % function [obj idx] = newItem(C,id,span[,{uicontrol properties}])
            % function [obj idx] = newItem(C,id,span,'panel')
            if nargin<4 || iscell(controlprop)
                if nargin<4, controlprop = {}; end
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
            idx = fn_find(id,{C.items.id});
            deleteValid([C.items(idx).obj])
            C.items(idx) = [];
            itemPositions(C)
        end
    end
    
    % Data (edit headers)
    methods
        function dataContextMenu(C)
            % init context menu
            m = C.contextmenu;
            delete(get(m,'children'))
            
            % create entries
            uimenu(m,'label','Edit header information', ...
                'callback',@(u,e)C.editHeader())
            uimenu(m,'label','Open data in a new xplor window', ...
                'callback',@(u,e)xplor(C.V.data))
            
            % make menu visible
            p = get(C.V.hf,'currentpoint'); p = p(1,1:2);
            set(m,'Position',p,'Visible','on')
        end
        function editHeader(C)
            data = C.V.data;
            curhead = data.header;
            newhead = xplr.editHeader(C.V.data);
            if isempty(newhead), return, end % user closed window: cancel
            dimchg = false(1,data.nd);
            for i=1:data.nd, dimchg(i) = ~isequal(newhead(i),curhead(i)); end
            if any(dimchg)
                dim = find(dimchg);
                newheaddim = xplr.dimheader(newhead(dim));
                C.V.data.updateData('chgdim',dim,[],data.data,newheaddim)
            end
        end
        function datachange(C,e)
            switch e.flag
                case 'global'
                    error 'global data change not handled'
                case 'chgdim'
                    % update dimension list
                    set(C.dimlist,'string',{C.V.data.header.label})
                otherwise
                    % no change needed
            end
        end
    end
    
    % Dimensions menu
    methods
        function dimensionContextMenu(C)
            % init context menu
            m = C.contextmenu;
            delete(get(m,'children'))
            
            % selected dimension(s)
            dim = get(C.dimlist,'value');
            dimID = [C.V.data.header(dim).dimID];
            
            % some strings to handle singular vs. plural
            dimstr = fn_switch(isscalar(dimID),'this dimension','these dimensions');
            filterstr = fn_switch(isscalar(dimID),'filter','filters');
            
            % add or change filter(s)
            % (2D with key 1)
            if length(dimID)==2
                % (using key 1)
                uimenu(m,'label','Filter with shared 2D filter','separator','on', ...
                    'callback',@(u,e)dimaction(C,'addfilter',{dimID},1))
                nextseparator = 'off';
            else
                nextseparator = 'off';
            end
            % (1D with key 1)
            uimenu(m, ...
                'label',['Filter with shared 1D ' filterstr],'separator',nextseparator, ...
                'callback',@(u,e)dimaction(C,'addfilter',num2cell(dimID),1))            
            % (more options: select among available keys)
            availablekeys = xplr.bank.availableFilterKeys('filterAndPoint');
            newkey = max(availablekeys)+1;
            keyvalues = [setdiff(availablekeys,1) newkey];
            m2 = uimenu(m,'label','Filter with');
            if length(dimID)==2
                for keyvalue=keyvalues
                    uimenu(m2,'label',['shared 2D filter ' num2str(keyvalue)], ...
                        'callback',@(u,e)dimaction(C,'addfilter',{dimID},keyvalue));
                end
            end
            for keyvalue=keyvalues
                uimenu(m2,'label',['shared 1D ' filterstr ' ' num2str(keyvalue)], ...
                    'callback',@(u,e)dimaction(C,'addfilter',num2cell(dimID),keyvalue));
            end
            uimenu(m2,'label',['private 1D ' filterstr], ...
                'callback',@(u,e)dimaction(C,'addfilter',num2cell(dimID),0))
            
            % remove filters in these dimensions
            uimenu(m,'label',['Remove ' filterstr],'separator','on', ...
                'callback',@(u,e)dimaction(C,'rmfilter',dimID))
            
            % filter all others dimension
            uimenu(m, ...
                'label',['View ' dimstr ', filter others'], 'separator', 'on', ...
                'callback',@(u,e)dimaction(C,'viewdim',dimID,1))
            uimenu(m,'label',['View ' dimstr ' in a new window'], ...
                'callback',@(u,e)dimaction(C,'newwindow_viewdim',dimID,1))
            
            % make menu visible
            p = get(C.V.hf,'currentpoint'); p = p(1,1:2);
            set(m,'Position',p,'Visible','on')
        end
        function dimaction(C,flag,dimID,varargin)
            % function dimaction(C,'addfilter',dimIDs[,key[,active]])
            % function dimaction(C,'rmfilter|showfilter',dimID)
            % function dimaction(C,'setactive',dimID,value)
            %---
            % if flag is 'addfilter', dims can be a cell array, to defined
            % several filters at once for example
            % dimaction(C,'addfilter',{[1 2] 3}) will add two filters,
            % first a 2D filter in dimensions [1 2], second a 1D filter in
            % dimension 3
            %
            % dimID is supposed to be the unique identifier of some
            % dimension(s), but for commodity it can also be the dimension
            % number, or the dimension label
            
            % other window
            if strfind(flag, 'newwindow') %#ok<STRIFCND>
                % open data in a new window: flag can be either
                % 'otherwindow' or 'otherwindow_action' where 'action' is
                % to be executed in this window
                V2 = xplor(C.V.data);
                tokens = regexp(flag, 'newwindow_(.*)','tokens');
                if ~isempty(tokens)
                    V2.C.dimaction(tokens{1}{1},dimID,varargin{:})
                end
                return
            end
            
            % commodity: convert dimension numbers or labels to dimension
            % identifiers
            dimID = C.V.data.dimensionID(dimID);
            
            % 'addfilter' flag -> several filters at once
            if strcmp(flag,'addfilter')
                % dims will be a cell array: list of dimensions, per filter
                % dimID will be an array: list of all affected dimensions
                if ~iscell(dimID)
                    if ~isscalar(dimID), error 'array of dimID values is ambiguous, use a cell array instead', end
                    dimIDs = {dimID};
                else
                    dimIDs = dimID; % several set of one or several dimensions
                    dimID = unique([dimIDs{:}]);
                end
                if length(dimID) < length([dimIDs{:}])
                    error 'some dimension is repeated in filter(s) definition'
                end
            end
            
            % list of filters in the selected dimensions
            filtersidx = find(fn_map({C.V.slicer.filters.dimID},@(dd)any(ismember(dd,dimID)),'array'));
            currentfiltersdim = C.V.slicer.filters(filtersidx); % current filters acting on dimensions within dd

            % filters to remove
            if ismember(flag,{'addfilter' 'rmfilter' 'viewdim'})
                % remove filter from the viewcontrol and the bank
                for filter = currentfiltersdim
                    C.removefilter(filter);
                end
                
                % remove filters from the slicer
                doslicing = strcmp(flag,'rmfilter'); % no need to reslice yet for 'addfilter', reslice will occur when adding the new filter(s)
                C.V.slicer.rmFilter(filtersidx, doslicing);
            end
            
            % filters to add
            if ismember(flag,{'addfilter' 'viewdim'})
                if nargin>=4, key = varargin{1}; else, key = 1; end
                if nargin>=5, active = varargin{2}; else, active = true; end
                if strcmp(flag,'addfilter')
                    adddimIDs = dimIDs; % already a cell array
                else
                    % add 1D filters il all dimensions that we do not want
                    % to view and that are not already filtered
                    noviewdimID = setdiff([C.V.data.header.dimID], dimID);
                    curfiltdimID = [C.V.slicer.filters.dimID];
                    adddimIDs = num2cell(setdiff(noviewdimID, curfiltdimID));
                end
                nadd = length(adddimIDs);
                if nadd > 0
                    if nadd>1 && isscalar(key), key = repmat(key,1,nadd); end
                    if nadd>1 && isscalar(active), active = repmat(active,1,nadd); end
                    % loop on dimension sets
                    newfilters = struct('dimID',cell(1,0),'F',[],'active',[]);
                    for i = 1:length(adddimIDs)
                        F = C.createFilterAndItem(adddimIDs{i},key(i),active(i));
                        newfilters(end+1) = struct('dimID',adddimIDs{i},'F',F,'active',active(i)); %#ok<AGROW>
                    end
                    C.V.slicer.addFilter({newfilters.dimID},[newfilters.F],[newfilters.active]) % slicing will occur now
                end
            end
            
            % show filter, set filter active
            switch flag
                case 'setactive'
                    active = varargin{1};
                    % show label(s) as enabled/disabled
                    for filter = currentfiltersdim
                        itemidx = fn_find({'filter' filter.dimID},{C.items.id});
                        hlab = C.items(itemidx).label;
                        set(hlab,'enable',fn_switch(active,'inactive','off'))
                        drawnow
                    end
                    % toggle filter active in slicer
                    C.V.slicer.chgFilterActive(filtersidx,active)
                case 'showfilter'
                    for filter = currentfiltersdim
                        F = filter.obj;
                        if ~isscalar(filter.dimID)
                            disp('cannot display list for ND filter')
                        elseif F.linkkey == 0
                            % private filter
                            combo = C.getPrivateLists();
                            combo.showList(F)
                        else
                            xplr.bank.showList(F);
                        end
                    end
            end

            % Empty the dimension selection
            set(C.dimlist,'value',[])
        end
        function addFilterItem(C,dimID,label,F,active)            
            % panel
            id = {'filter' dimID};
            [panel, itemidx] = C.newItem(id,1,'panel');
            backgroundColor = xplr.colors('linkkey',F.linkkey);
            panel.BackgroundColor = backgroundColor;
            
            % store the filter
            C.items(itemidx).F = F;
            
            % label
            hlab = uicontrol('parent',panel, ...
                'position',[20 5 300 15], ...
                'style','text','string',label,'horizontalalignment','left', ...
                'backgroundcolor',backgroundColor, ...
                'enable', fn_switch(active,'inactive','off'), ...
                'buttondownfcn',@(u,e)clickFilterItem(C,dimID,id));
            C.items(itemidx).label = hlab;
            
            % buttons
            [ii, jj] = ndgrid(-2:2);
            x = min(1,abs(abs(ii)-abs(jj))*.5);
            x(x==1) = NaN; x = repmat(x,[1 1 3]);
            
            % cross button to remove the filter
            rmFilterButton = uicontrol('parent',panel,'cdata',x, ...
                'unit', 'normalized', ...
                'position', [ 0.95 0.5 0.05 0.5 ], ...
                'callback',@(u,e)C.dimaction('rmfilter',dimID));
            fn_controlpositions(rmFilterButton, panel, [1 .5 0 .5], [-11 0 11 0]);
            
            % checkbox to disable and enable the filter
            uicontrol('parent',panel, ...
                'backgroundcolor',backgroundColor, ...
                'Style','checkbox', 'Value',active, ...
                'position', [ 6 6 13 12 ], ...
                'callback',@(u,e)C.dimaction('setactive',dimID,get(u,'value')));
            
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
            idxitem = fn_find(id, {C.items.id});
            idx0 = idxitem-(idxfilter(1)-1);
            idxother = setdiff(1:nfilter,idx0);
            obj = C.items(idxitem).obj;
            pos0 = get(obj,'position');
            ystep = 24;
            
            % move
            p0 = get(hf,'currentpoint'); p0 = p0(1,2); % only vertical position matters
            newidx = [];
            moved = fn_buttonmotion(@move,hf,'moved?');
            function move
                p = get(hf,'currentpoint'); p = p(1,2);
                newidx = fn_coerce( idx0 - round((p-p0)/ystep), 1, nfilter);
                % set all items position
                C.items(idxfilter) = filteritems([idxother(1:newidx-1) idx0 idxother(newidx:end)]);
                C.itemPositions
                % set selected item position
                newpos = pos0; newpos(2) = pos0(2) + fn_coerce(p-p0,[idx0-nfilter idx0-1]*ystep);
                set(obj,'position',newpos)
            end
            if moved
                % re-position correctly the selected item
                C.itemPositions
                % apply filters permutation
                perm = [idxother(1:newidx-1) idx0 idxother(newidx:end)];
                C.V.slicer.permFilters(perm)
            end
            
            % show filter if there was no move
            if ~moved, dimaction(C,'showfilter',d), end
        end
    end
    
    % Filters display
    methods (Access='private')
        function removefilter(C,filter)
            % remove the filter from the viewcontrol and the bank
            % this function does not remove the filter from the slicer
            
            % if filter is empty, does nothing and leave the function
            if isempty(filter), return, end
            % remove filter from the items
            C.rmItem({'filter' filter.dimID})
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
        function F = createFilterAndItem(C,dimID,key,active)
            % create filter or get existing one from the
            % related public filters set
            
            header = C.V.data.headerByID(dimID);
            % if the filter has to be private
            if key == 0
                % create private filter
                F = xplr.filterAndPoint(header);
                % show filter in combo
                if isscalar(dimID)
                    combo = C.getPrivateLists();
                    if active, combo.showList(F), end
                end
            else
                % search for the filter in the bank with key and dimension
                doshow = true;
                F = xplr.bank.getFilterAndPoint(key,header,C,doshow);
            end
            
            % add the filter to the items, it is important that
            % filter.linkkey is set before using addFilterItem
            % TODO: change how the string filter is shifted
            str = ['filter ' header.label ' (' char(F.F.slicefun) ')'];
            C.addFilterItem(dimID,str,F,active)
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
                combo = xplr.listcombo(C.V.panels.listcombo);
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
    
    
end