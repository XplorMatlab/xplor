classdef list < xplr.graphnode
    % function L = list(filter[,'in',uipanel][,other options...])
        
    properties (SetAccess='private')
        F           % xplr.filterAndPoint object
        seltype
        hlist
        hp  % direct parent
        hf  % figure parent
        hlabel
        menu
        valuestr    % precomputed list of values
    end
    properties (SetObservable, AbortSet=true)
        selmultin = true;
        scrollwheel = 'on'; % 'on', 'off' or 'default'
        selectionpromptname = 'none'; % 'all', 'groups' or 'none'
    end
    
    % Constructor and Destructor
    methods
        function L = list(F,varargin)
            % options for initialization
            opt = struct( ...
                'in',                   [] ...
                );
            if nargin==0
                head = xplr.header('testheader',10);
                F = xplr.filterAndPoint(head);
            end
            [opt, optadd] = parseInput(opt,varargin{:});
            
            % check filter
            L.F = F;
            if F.ndin~=1 || F.ndout~=1, error 'input and output of list filter must be one-dimensional', end
            if ~isa(F,'xplr.filterAndPoint'), error 'list can act only on a filterAndPoint object', end
            L.seltype = fn_switch(F.headerin.categorical,'indices','point1D');
            
            % watch filter deletion
            addlistener(L.F,'ObjectBeingDestroyed',@(u,e)delete(L));
            
            % add 'SoftSelection' label to output header
            F.augmentHeader('SoftSelection','logical')
            
            % uipanel container
            newfigure = isempty(opt.in);
            if newfigure
                % create uipanel in a new figure
                L.hp = uipanel('parent',figure);
            elseif strcmp(get(opt.in,'type'),'uipanel')
                L.hp = opt.in;
            else
                error 'input container must be an uipanel object'
            end
            
            % parent figure
            L.hf = fn_parentfigure(L.hp);

            % create several components
            % (list)
            L.hlist = uicontrol('parent',L.hp,'style','listbox','min',0,'max',2, ...
                'callback',@(hlist,evnt)event(L,'select'), ...
            	'keypressfcn',@(hlist,evnt)keypress(L,evnt));
            fn_controlpositions(L.hlist,L.hp,[0 0 1 1],[8 5 -16 -5-21-2])
            % (label)
            L.hlabel = uicontrol('parent',L.hp,'style','text', ...
                'string',L.F.headerout.label, ...
                'horizontalalignment','center', ...
                'backgroundcolor',xplr.colors('linkkey',L.F.linkkey));
            fn_controlpositions(L.hlabel,L.hp,[0 1 1 0],[8 -21 -8-18 18])
            % (close button)
            if ~newfigure
                x = fn_printnumber(ones(18),'x','pos','center')';
                x(x==1) = NaN; x = repmat(x,[1 1 3]);
                hclose = uicontrol('parent',L.hp,'cdata',x,'callback',@(u,e)delete(L));
                fn_controlpositions(hclose,L.hp,[1 1],[-8-18 -3-18 18 18])
            end
            % (group button)
            ctrl = fn_propcontrol(L,'selmultin','togglebutton', ...
                {'parent',L.hp,'string','G'});
            fn_controlpositions(ctrl.hu,L.hp,[1 1],[-8-38 -3-18 18 18])
            
            % context menu
            initlocalmenu(L)
            
            % list and event (bottom-up)
            if boolean(L.scrollwheel)
                L.scrollwheel = 'on'; % this will automaticall register scroll wheel
            end
            if isempty(get(L.hf,'WindowButtonMotionFcn'))
                % force update of current position when moving the mouse
                % around
                set(L.hf,'WindowButtonMotionFcn',@(u,e)donothing())
            end
            
            % watch filter
            function filterchanged(~,e)
                if strcmp(e.type,'filter'), displayselection(L), end
            end
            connectlistener(F,L,'ChangedOperation',@filterchanged);
            
            % auto-delete
            set(L.hlist,'deletefcn',@(u,e)delete(L))
            addlistener(F,'ObjectBeingDestroyed',@(u,e)delete(L));

            % update display (here, just sets the correct value)
            preformatvalues(L)
            displayselection(L)
            
            % set more properties
            if ~isempty(optadd)
                set(L,optadd{:})
            end
            
            % put object in base workspace for debugging purposes
            if nargin==0, assignin('base','L',L), end
        end
        function initlocalmenu(L)
            delete(L.menu)
            L.menu = uicontextmenu('parent',L.hf);
            m = L.menu;
            set(L.hlist,'UIContextMenu',m)
            
            uimenu(m,'label','New singleton selections','callback',@(u,e)event(L,'newuni'))
            uimenu(m,'label','New group selection [A]','callback',@(u,e)event(L,'newgroup'))
            uimenu(m,'label','Add to selection','callback',@(u,e)event(L,'add'))

            uimenu(m,'label','Define new group...','callback',@(u,e)event(L,'definegroup'),'separator','on')

            uimenu(m,'label','Sort selections according to list order','callback',@(u,e)event(L,'sort'),'separator','on')
            uimenu(m,'label','Reorder selections...','callback',@(u,e)event(L,'reorder'))
            
            uimenu(m,'label','Remove highlighted group(s)','callback',@(u,e)event(L,'rmgroup'),'separator','on')
            uimenu(m,'label','Remove all groups','callback',@(u,e)event(L,'rmgroupall'))
            uimenu(m,'label','Remove highlighted individuals','callback',@(u,e)event(L,'rmuni'))
            uimenu(m,'label','Remove all individuals','callback',@(u,e)event(L,'rmuniall'))
            uimenu(m,'label','Remove highlighted selections','callback',@(u,e)event(L,'rm'))
            uimenu(m,'label','Remove all selections','callback',@(u,e)event(L,'rmall'))

            uimenu(m,'label','Select all','separator','on','callback',@(u,e)event(L,'selectall'))
            
            fn_propcontrol(L,'selmultin', ...
                {'menuval', {true false}, {'individuals' 'group'}}, ...
                {'parent',m,'label','Temporary selection','separator','on'});
            fn_propcontrol(L,'selectionpromptname', ...
                {'menu' {'all' 'groups' 'none'} {'all selections' 'group selections only' 'none'}}, ...
                {'parent',m,'label','Prompt for selection name'});
            
            % scroll wheel behavior: changing it would make sense only if
            % list is inside a figure with other elements, which does not
            % occur in XPLOR at the present time
            % and anyway, there are some bugs, and the whole
            % windowcallbackmanager thing needs to be replaced by calls to
            % iptaddcallback
            %             m1 = uimenu(m,'label','scroll wheel','separator','on');
            %             fn_propcontrol(L,'scrollwheel', ...
            %                 {'menu', {'on' 'off' 'default'}}, ...
            %                 {'parent',m1,'label','Scroll wheel behavior'});
            %             uimenu(m1,'label','make default in figure', ...
            %                 'callback',@(u,e)set(L,'scrollwheel','default'));
        end
        function delete(L)
            delete@xplr.graphnode(L)
            if ~isvalid(L) && ~isprop(L,'hlist'), return, end
            deleteValid(L.hlist,L.hlabel)
        end
    end
       
    % Events
    methods
        function keypress(L,e)
            switch e.Key
                case 'a'
                    event(L,'newgroup')
                case {'insert' 'delete'}
                    event(L,'scroll',fn_switch(e.Key,'insert',-1,'delete',1))
            end
        end
        function event(L,flag,varargin)
            % possible values for flag:
            % - select          selection by user click in list display
            % - unisel          selection by user double-click
            % - ('scroll',n)    scrolling
            % - newuni, newgroup
            % - definegroup
            % - add
            % - sort
            % - reorder...
            % - rmall, rmgroup, rmgroupall, rmuni, rmuniall, rm         
            
            % selected list entries
            val = get(L.hlist,'value');
            %L.F.shared.list.cursel = val; % make available to all lists acting on this filter what is the new entries selection
            
            % get the current selections
            selinds = L.F.F.indices; % it is important here not to get L.F.indices, because we do not want the point selection to appear here
            nsel = length(selinds);
            isunisel = false(1,nsel); 
            for i=1:nsel, isunisel(i) = isscalar(selinds{i}); end % faster than calling fn_map
            
            % which selections are soft
            softsel = L.F.F.headerout.getValue('SoftSelection'); % same as above
            softsel = [softsel{:}]; % faster than cell2mat
            if isempty(softsel), softsel = false(1,nsel); end
            isoft = find(softsel);
            nsoft = length(isoft);
            isolid = find(~softsel);
            nsolid = nsel-nsoft;
            
            % modify flag if needed
            switch flag
                case 'select'
                    if strcmp(get(L.hf,'selectiontype'),'open')
                        flag = 'unisel';
                    end
                case 'selectall'
                    set(L.hlist,'value',1:L.F.szin)
                    flag = 'select';
                case 'scroll'
                    n = varargin{1};
                    soft = selinds(softsel);
                    if length(soft)>=2 && all(fn_map(@isscalar,soft)) ...
                            && all(ismember(diff([soft{:}]),[0 1]))
                        % current soft selection consists of a range of
                        % consecutive values -> select the same number of
                        % values before or after
                        % Note that when hitting the lower or upper value,
                        % this value can be repeated several times in the
                        % selection in order to remember how many values
                        % where selected in total.
                        soft = [soft{:}];
                        nval = length(soft);
                        if n<0
                            % decreasing values
                            if any(diff(soft)==0)
                                % if last value was selected several times
                                % show the nval last values
                                val = soft(end) + (-nval+1:0);
                            else
                                % otherwise show the nval values before
                                % soft(1)
                                val = soft(1) + n*nval + (0:nval-1);
                            end
                        else
                            % increasing values, same idea
                            if any(diff(soft)==0)
                                val = soft(1) + (0:nval-1);
                            else
                                val = soft(end) + (n-1)*nval + (1:nval);
                            end
                        end
                        val = fn_coerce(val,1,L.F.szin);
                    else
                        val = fn_coerce(L.F.index+n,1,L.F.szin);
                        if val==L.F.index, return, end
                    end
                    set(L.hlist,'value',val,'listboxtop',val(1)-2);
                    flag = 'select';
            end            
            
            % action (or only determine shich selections to remove and
            % which to add)
            idxrm = []; newsel = []; newissoft = false;
            switch flag
                case 'select'
                    if isscalar(val)
                        L.F.index = val;
                    end
                    if isempty(val) || (isscalar(val) && isempty(isolid)) ...
                            || (isscalar(val) && isequal({val},selinds(isoft)))
                        % remove temporary selection in the following
                        % cases:
                        % - user unselected all list items
                        % - no solid selection and a single selected item
                        % - repeated selection of temporaray selection item
                        idxrm = isoft;
                    else
                        % new temporary selection
                        idxrm = isoft;
                        newsel = buildCurrentSelection(L,L.selmultin);
                        newissoft = true;
                    end
                case 'unisel'
                    % double-click -> make new solid selection with current
                    % index, or remove it
                    if ~isscalar(val), return, end
                    kunisel = find(isunisel);
                    idxrm = kunisel([selinds{kunisel}]==val); % index of already-existing selection with this value
                    if isempty(idxrm)
                        % create
                        newsel = buildCurrentSelection(L,true);
                    elseif softsel(idxrm)
                        % make solid
                        newsel = buildCurrentSelection(L,true);
                    else
                        % remove
                    end
                case 'newuni'
                    idxrm = isoft;
                    newsel = buildCurrentSelection(L,true);
                case 'newgroup'
                    idxrm = isoft;
                    newsel = buildCurrentSelection(L,false);
                case 'definegroup'
                    str = inputdlg('Define selection','',1,{['1:' num2str(L.F.szin)]}); % TODO: continue!!!
                    if isempty(str), disp 'interrupted', return, end
                    try
                        val = evalin('base',['[' str{1} ']']);
                        if ~iscell(val), val = {val}; end
                        newsel = xplr.selectionnd(length(val));
                        for i=1:length(val)
                            newsel(i) = xplr.selectionnd(L.seltype,val{i},L.F.headerin.n);
                        end
                    catch
                        errordlg('Command could not be evaluated correctly')
                        return
                    end
                case 'add'
                    newsel = xplr.selectionnd('point1D',val);
                    if nsoft, updateSelection(L.F,'remove',isoft), end % remove all soft selections
                    if ~isempty(isolid)
                        updateSelection(L.F,'add',isolid(end),newsel)
                    else
                        updateSelection(L.F,'new',newsel)
                    end
                    return
                case 'sort'
                    % remove all soft selections
                    if nsoft, updateSelection(L.F,'remove',isoft), end
                    selinds(isoft) = [];
                    % reorder other selections
                    idxfirst = fn_map(selinds,@(v)v(1),'array');
                    [~, ord] = sort(idxfirst);
                    updateSelection(L.F,'perm',ord)
                    return
                case 'reorder'
                    % first remove all soft selections
                    if nsoft, updateSelection(L.F,'remove',isoft), end
                    nsel = nsolid;
                    % prompt for reordering
                    ord = fn_input('new order',1:nsel);
                    if length(unique(ord))<length(ord) || ~all(ismember(ord,1:nsel))
                        waitfor(errordlg('Not a valid permutation or subset'))
                        return
                    end
                    if length(ord)<nsel
                        answer = questdlg('This is not a permutation: some selections will be removed','', ...
                            'OK','Cancel','OK');
                        if strcmp(answer,'Cancel'), return, end
                        idxrm = setdiff(1:nsel,ord);
                        updateSelection(L.F,'remove',idxrm)
                        [~, ordsort] = sort(ord);
                        ord(ordsort) = 1:length(ord);
                    end
                    updateSelection(L.F,'perm',ord)
                    return
                case 'rmall'
                    % set empty selections rather than remove all existing
                    % ones: this can performs some clean-up when errors
                    % occured previously
                    newsel = xplr.selectionnd.empty(1,0);
                    updateSelection(L.F,'all',newsel)
                    return
                case {'rmgroup' 'rmgroupall' 'rmuni' 'rmuniall' 'rm'}
                    if strfind(flag,'all')
                        range = 1:nsel;
                        flag = strrep(flag,'all','');
                    else
                        range = fn_find(@(x)intersect(x,val),selinds);
                    end
                    rmmask = false(1,nsel);
                    switch flag
                        case 'rmgroup'
                            rmmask(range) = ~isunisel(range);
                        case 'rmuni'
                            rmmask(range) = isunisel(range);
                        case 'rm'
                            rmmask(range) = true;
                    end
                    rmmask(isoft) = true; % in any case, remove all soft selections
                    idxrm = find(rmmask);
            end
            
            % prompt for name of new selections
            name_options = {};
            if ~strcmp(L.selectionpromptname,'none') && ~newissoft
                anyname = false;
                nnew = length(newsel);
                names = cell(1,nnew);
                for i = 1:nnew
                    if strcmp(L.selectionpromptname,'groups') && isscalar(newsel(i).dataind), continue, end
                    name = inputdlg('Group name','xplor');
                    if isempty(name), continue, end
                    names{i} = name;
                    anyname = true;
                end
                if anyname
                    name_options = {'Name' names};
                end
            end
            
            % remove/change/add new selections
            nnew = length(newsel);
            nrm = length(idxrm);
            if nnew==0 && nrm==0
                % happens for example when there are no selection, and the
                % point is moved -> nothing to do
            elseif nrm==0
                updateSelection(L.F,'new',newsel,'SoftSelection',newissoft)
            elseif nnew==0
                updateSelection(L.F,'remove',idxrm)
            elseif nnew==nrm
                updateSelection(L.F,'chg',idxrm,newsel,'SoftSelection',newissoft,name_options{:})
            elseif nnew>nrm
                updateSelection(L.F,'chg&new',{idxrm nsel+(1:nnew-nrm)},newsel,'SoftSelection',newissoft,name_options{:})
            elseif nrm>nnew
                updateSelection(L.F,'chg&rm',{idxrm(1:nnew) idxrm(nnew+1:nrm)},newsel,'SoftSelection',newissoft,name_options{:})
            end
        end
        function sel = buildCurrentSelection(L,domultin)
            val = get(L.hlist,'value');
            if isempty(val), sel = []; return, end
            if domultin 
                nsel = length(val);
                sel = xplr.selectionnd(nsel);
                for i=1:nsel
                    sel(i) = xplr.selectionnd(L.seltype,val(i),L.F.headerin.n);
                end
            elseif L.F.headerin.ismeasure && ~isscalar(val) && all(diff(val)==1)
                % selection is a segment rather than a mere list of points
                sel = xplr.selectionnd('line1D',val([1 end]) + [-.5 .5],L.F.headerin.n);
            else
                sel = xplr.selectionnd(L.seltype,val,L.F.headerin.n);
            end
        end
    end
        
    % Get/Set - scroll wheel
    methods
        function set.scrollwheel(L,flag)
            switch flag
                case 'on'
                    fn_scrollwheelregister(L.hlist,@(n)event(L,'scroll',n)) %#ok<MCSUP>
                    L.scrollwheel = 'on';
                case 'default'
                    fn_scrollwheelregister(L.hlist,@(n)event(L,'scroll',n),'default') %#ok<MCSUP>
                    L.scrollwheel = 'on';
                case 'off'
                    fn_scrollwheelregister(L.hlist,flag) %#ok<MCSUP>
                    L.scrollwheel = 'off';
                otherwise
                    error 'scrollwheel value must be ''off'', ''on'' or ''default'''
            end
        end
    end
    
    % Get/Set
    methods
        function set.selmultin(L,val)
            if val==L.selmultin, return, end
            L.selmultin = val;
            % update selection
            event(L,'select')
        end
    end
    
    % Display
    methods (Access = 'private')
        function preformatvalues(L)
            L.valuestr = L.F.headerin.getItemNames();
        end
        function displayselection(L)
            % init list with names of items
            str = L.valuestr;
            
            selinds = L.F.F.indices; % it is important here not to get L.F.indices, because we do not want the point selection to appear here
            nsel = length(selinds);
            isunisel = false(1,nsel); 
            for i=1:nsel, isunisel(i) = isscalar(selinds{i}); end % faster than calling fn_map
            softsel = L.F.F.headerout.getValue('SoftSelection'); % same as above
            softsel = [softsel{:}]; % faster than cell2mat
            if isempty(softsel), softsel = false(1,nsel); end
                        
            % mark selections
            for ksel=find(isunisel & ~softsel)
                ind = selinds{ksel};
                str{ind} = [str{ind} '[' num2str(ksel) ']'];
            end
            for ksel=find(~isunisel & ~softsel)
                for ind = selinds{ksel}
                    str{ind} = [str{ind} '[group' num2str(ksel) ']'];
                end
            end
            for ksel=find(isunisel & softsel)
                ind = selinds{ksel};
                str{ind} = [str{ind} '*'];
            end
            for ksel=find(~isunisel & softsel)
                for ind = selinds{ksel}
                    str{ind} = [str{ind} '*g'];
                end
            end
            
            % update display!
            top = get(L.hlist,'listboxtop'); % the portion of the list that is shown moves when setting the string -> reset it to its current position
            set(L.hlist,'string',str,'ListboxTop',top)
            if nsel==0
                % if no selection, highlight the point selection
                set(L.hlist,'value',L.F.index)
            else
                set(L.hlist,'value',[selinds{softsel}])
            %             elseif isfield(L.F.shared,'list')
            %                 set(L.hlist,'value',L.F.shared.list.cursel)
            end
        end
            
    end
    
end
     
%---
function donothing()
% setting this function as main figure WindowButtonMotionFcn forces the
% figure to update CurrentPoint whenever the mouse is moved, and therefore
% have this property set currently even when clicking active controls
% (unfortunately, clicking a control does not set this property)
end
   