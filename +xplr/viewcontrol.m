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
            datastr = 'Data';
            if ~isempty(V.data.name), datastr = [datastr ' (' V.data.name ')']; end
            C.new_item('data',1, ...
                {'style','text','string',V.data.name, ...
                'backgroundcolor',xplr.colors('gui.controls.dataname'), ...
                'enable','inactive','buttondownfcn',@(u,e)C.dataContextMenu()})
            % (list of data dimensions)
            C.dimlist = C.new_item('dimlist',4, ...
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
            fn_pixelsizelistener(C.hp,@(u,e)item_positions(C))
            
            % note that other fields will be added, e.g. in addFilterItem
            C.items = struct('id',cell(1,0),'span',[],'obj',[]);
        end
        function item_positions(C,idx)
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
        function [obj, idx] = new_item(C,id,span,controlprop)
            % function [obj idx] = new_item(C,id,span[,{uicontrol properties}])
            % function [obj idx] = new_item(C,id,span,'panel')
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
            item_positions(C,idx)
            if nargout==0, clear obj, end
        end
        function item = get_item(C,id)
            idx = fn_find(id,{C.items.id});
            item = C.items(idx);
        end
        function remove_item(C,id)
            idx = fn_find(id,{C.items.id});
            deleteValid([C.items(idx).obj])
            C.items(idx) = [];
            item_positions(C)
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
    
    % Dimensions menu and actions
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
            
            % convert dimension numbers or labels to dimension identifiers
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
                    C.remove_filter_item(filter.dimID);
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
                    dimIDs_add = dimIDs; % already a cell array
                else
                    % add 1D filters il all dimensions that we do not want
                    % to view and that are not already filtered
                    noviewdimID = setdiff([C.V.data.header.dimID], dimID, 'stable');
                    curfiltdimID = [C.V.slicer.filters.dimID];
                    dimIDs_add = setdiff(noviewdimID, curfiltdimID, 'stable');
                    % among these dimensions, attempt to find pairs of
                    % measure headers with same units to set 2D filter
                    % instead of two 1D filters
                    head = C.V.data.headerByID(dimIDs_add);
                    connections = measure_grouping(head);
                    pairs = {};
                    while any(connections(:))
                        [i, j] = find(connections,1,'first');
                        pairs{end+1} = dimIDs_add(sort([i j]));
                        connections([i j],:) = false;
                        connections(:, [i j]) = false;
                    end
                    dimIDs_add = [pairs num2cell(setdiff(dimIDs_add,[pairs{:}],'stable'))];
                end
                nadd = length(dimIDs_add);
                if nadd > 0
                    if nadd>1 && isscalar(key), key = repmat(key,1,nadd); end
                    if nadd>1 && isscalar(active), active = repmat(active,1,nadd); end
                    % loop on dimension sets
                    newfilters = struct('dimID',cell(1,0),'F',[],'active',[]);
                    for i = 1:length(dimIDs_add)
                        F = C.create_filter_and_item(dimIDs_add{i},key(i),active(i));
                        newfilters(end+1) = struct('dimID',dimIDs_add{i},'F',F,'active',active(i)); %#ok<AGROW>
                    end
                    C.V.slicer.addFilter({newfilters.dimID},[newfilters.F],[newfilters.active]) % slicing will occur now
                else
                    % we might have removed filters before without updating
                    % completely the slice
                    C.V.slicer.applyPending()
                end
                
                % adjust display mode and layout if it seems appropriate
                D = C.V.D;
                if strcmp(flag,'viewdim')
                    if isscalar(dimID)
                        D.set_dim_location(dimID,'x',strcmp(D.displaymode,'time courses'))
                        D.displaymode = 'time courses';
                    elseif length(dimID)==2
                        D.set_dim_location(dimID,{'x' 'y'},strcmp(D.displaymode,'image'))
                        D.displaymode = 'image';
                    end
                else
                    nsdimID = non_singleton_dimID(C.V.slice.header);
                    if isscalar(nsdimID)
                        D.set_dim_location(nsdimID,'x',strcmp(D.displaymode,'time courses'))
                        D.displaymode = 'time courses';
                    end
                end
            end
            
            % show filter, set filter active
            switch flag
                case 'setactive'
                    active = varargin{1};
                    % show label(s) as enabled/disabled
                    for filter = currentfiltersdim
                        item = C.get_item({'filter' filter.dimID});
                        hlab = [item.filter_label item.dimension_label];
                        set(hlab,'enable',fn_switch(active,'inactive','off'))
                        set(item.checkbox,'value',active)
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
    end
    
    % Filters display
    methods (Access='private')
        function remove_filter_item(C,dimID)
            % remove the filter from the viewcontrol and the bank
            % this function does not remove the filter from the slicer
            
            % if filter is empty, does nothing and leave the function
            if isempty(dimID), return, end
            % get filter
            id = {'filter' dimID};
            F = C.get_item(id).F;
            % remove filter from the items
            C.remove_item(id)
            % remove filter from the lists display
            % if the filter is private
            if F.linkkey == 0
                % remove the filter from the combo
                combo = C.getPrivateLists();
                combo.removeList(F)
            else
                % viewcontrol object C will be unregistered for the users
                % list of filter F; if this list will become empty, F will
                % be unregistered from the filters set
                xplr.bank.unregisterFilter(F,C) 
            end
        end
        function F = create_filter_and_item(C,dimID,key,active,show_new_filter)
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
                if nargin<5, show_new_filter = true; end
                F = xplr.bank.getFilterAndPoint(key,header,C,show_new_filter);
            end
            
            % panel
            id = {'filter' dimID};
            [panel, itemidx] = C.new_item(id,1,'panel');
            backgroundColor = xplr.colors('linkkey',F.linkkey);
            panel.BackgroundColor = backgroundColor;
            
            % store the filter
            C.items(itemidx).F = F;
            
            % filter and dimension labels
            % (create labels)
            filter_label_name = uicontrol('parent',panel, ...
                'style','text','string','filter','horizontalalignment','left', ...
                'backgroundcolor',backgroundColor, ...
                'enable', fn_switch(active,'inactive','off'), ...
                'buttondownfcn',@(u,e)click_filter_item(C,dimID), ...
                'uicontextmenu',uicontextmenu(C.V.hf,'callback',@(m,e)F.context_menu(m)));
            dimension_label = uicontrol('parent',C.hp, ...
                'style','text','horizontalalignment','left', ...
                'string',fn_strcat({F.headerin.label},'-'), ...
                'buttondownfcn',@(u,e)move_filtered_dimension(C,dimID), ...
                ... 'buttondownfcn',@(u,e)click_filter_item(C,dimID,id), ...
                'enable', fn_switch(active,'inactive','off'));
            filter_label_op = uicontrol('parent',panel, ...
                'style','text','string',['(' F.F.slicefunstr ')'],'horizontalalignment','left', ...
                'backgroundcolor',backgroundColor, ...
                'enable', fn_switch(active,'inactive','off'), ...
                'buttondownfcn',@(u,e)click_filter_item(C,dimID), ...
                'uicontextmenu',uicontextmenu(C.V.hf,'callback',@(m,e)F.context_menu(m)));
            % (adjust their positions based on their extents)
            w_name = filter_label_name.Extent(3);
            w_dim = dimension_label.Extent(3);
            set(filter_label_name,'position',[20 5 w_name 15])
            fn_controlpositions(dimension_label,panel,[],[20+w_name-1 5-1 w_dim 15])
            set(filter_label_op,'position',[20+w_name+w_dim 5 300 15])
            % (store handles)
            C.items(itemidx).filter_label = [filter_label_name filter_label_op];
            C.items(itemidx).dimension_label = dimension_label;
            
            % change filter label upon operation change
            function check_operation_change(~,e)
                if strcmp(e.type,'operation')
                    label = ['(' F.F.slicefunstr ')'];
                    set(filter_label_op,'string',label)
                end
            end
            hl = addlistener(F.F,'ChangedOperation',@check_operation_change);
            addlistener(filter_label_op,'ObjectBeingDestroyed',@(u,e)delete(hl));
            
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
            C.items(itemidx).checkbox = uicontrol('parent',panel, ...
                'backgroundcolor',backgroundColor, ...
                'Style','checkbox', 'Value',active, ...
                'position', [ 6 6 13 12 ], ...
                'callback',@(u,e)C.dimaction('setactive',dimID,get(u,'value')));
            
        end
        function click_filter_item(C,dimID)
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
            id = {'filter' dimID};
            idxitem = fn_find(id, {C.items.id});
            idx0 = idxitem-(idxfilter(1)-1);
            idxother = setdiff(1:nfilter,idx0);
            obj = C.items(idxitem).obj;
            pos0 = get(obj,'position');
            ystep = 24;
            
            % move
            p0 = get(hf,'currentpoint'); p0 = p0(1,2); % only vertical position matters
            newidx = [];
            moved = fn_buttonmotion(@move,hf,'moved?','pointer','hand');
            function move
                p = get(hf,'currentpoint'); p = p(1,2);
                newidx = fn_coerce( idx0 - round((p-p0)/ystep), 1, nfilter);
                % set all items position
                C.items(idxfilter) = filteritems([idxother(1:newidx-1) idx0 idxother(newidx:end)]);
                C.item_positions
                % set selected item position
                newpos = pos0; newpos(2) = pos0(2) + fn_coerce(p-p0,[idx0-nfilter idx0-1]*ystep);
                set(obj,'position',newpos)
            end
            if moved
                % re-position correctly the selected item
                C.item_positions
                % apply filters permutation
                perm = [idxother(1:newidx-1) idx0 idxother(newidx:end)];
                C.V.slicer.permFilters(perm)
            end
            
            % show filter if there was no move
            if ~moved, dimaction(C,'showfilter',dimID), end
        end
    end
    
    % Fancy moving dimensions from the filter items to the graph and
    % vice-versa
    methods
        % moving from the filters to the graph: when this happens, finishes
        % by a call to xplr.displaylabels.labelMove to select where to
        % locate the dimension in the graph
        function move_filtered_dimension(C,dimID)
            
            % move dimension label only if the filter is active
            id = {'filter' dimID};
            itemidx = fn_find(id,{C.items.id});
            item = C.items(itemidx);
            label = item.dimension_label;
            active = boolean(item.checkbox.Value);
            if ~active, return, end
            
            % move
            hf = fn_parentfigure(C.hp);
            p0 = get(hf,'currentpoint'); p0 = p0(1,1:2);
            pos0 = get(label,'pos');
            controls_width = C.hp.Position(3);
            panel = [];
            
            fn_buttonmotion(@movesub,hf,'pointer','hand')
            function movesub
                % once filter has been removed, do not execute this
                % callback any more
                if ~active, return, end
                
                % move label
                p = get(hf,'currentpoint'); p = p(1,1:2);
                pos = pos0; pos(1:2) = pos0(1:2)+(p-p0);
                set(label,'pos',pos)
                
                % disable filter if we exited the panel by the right side,
                % and immediately run displaylabel 'labelMove' method to
                % allow choosing where to position the dimension!!
                active = (p(1) <= controls_width);
                if ~active
                    % stop filtering
                    % prevent panel of being deleted now as this leads to a
                    % strange bug when label is deleted or even only hidden!!
                    panel = item.obj;
                    set(panel,'visible','off')
                    C.items(itemidx).obj = []; 
                    C.dimaction('rmfilter',dimID)
                    % move dimension label inside graph (note that this
                    % will call fn_buttonmotion in the same figure, and
                    % therefore terminate the current fn_buttonmotion)
                    % we activate immediate display update
                    if isscalar(dimID)
                        L = C.V.D.labels;
                        mem_do_update = L.doImmediateDisplay;
                        L.doImmediateDisplay = true;
                        L.labelMove(dimID,false)
                        L.doImmediateDisplay = mem_do_update;
                    end
                end
            end
            
            % and put back at original position when we release the mouse
            % button!
            if active
                set(label,'pos',pos0)
            else
                delete(label)
                delete(panel)
            end
        end
        % moving from the graph to the filters: the methods below will be
        % called by xplr.displaylabels.labelMove
        function show_inoperant_filter(C,dimID)
            create_filter_and_item(C,dimID,1,true,false);
        end
        function activate_inoperant_filter(C,dimID)
            item = C.get_item({'filter' dimID});
            C.V.slicer.addFilter(dimID,item.F) % slicing will occur
        end
        function remove_inoperant_filter(C,dimID)
            remove_filter_item(C,dimID)
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