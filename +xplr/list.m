classdef list < xplr.graphnode
    % function L = list(filter[,'in',fig/axes/uicontrol/ctrl+label][,other options...])
        
    properties (SetAccess='private')
        F           % xplr.filterAndPoint object
        doroi
        hu
        hp  % direct parent
        hf  % figure parent
        hlabel
        menu
        menuitems
        valuestr    % precomputed list of values
    end
    properties
        selmultin = true;
        scrollwheel = 'on'; % 'on', 'off' or 'default'
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
                F = xplr.filterAndPoint(head,'indices');
            end
            [opt optadd] = parseInput(opt,varargin{:});
            
            % check filter
            L.F = F;
            if F.ndin~=1 || F.ndout~=1, error 'input and output of list filter must be one-dimensional', end
            if ~isa(F,'xplr.filterAndPoint'), error 'list can act only on a filterAndPoint object', end
            L.doroi = strcmp(F.F.type,'selection');
            
            % add 'SoftSelection' label to output header
            F.augmentHeader('SoftSelection','logical')
            
            % figure and axes
            if isempty(opt.in), opt.in = gcf; end
            if length(opt.in)==2
                L.hlabel = opt.in(2); % handle for label is given in input; label will be properly displayed later
                opt.in = opt.in(1);
            end
            if ~ishandle(opt.in) && mod(opt.in,1)==0 && opt.in>0, figure(opt.in), end
            switch get(opt.in,'type')
                case 'figure'
                    L.hp = opt.in;
                    figure(opt.in), set(L.hp,'menubar','none'), delete(findall(opt.in,'parent',opt.in))
                    L.hu = uicontrol('units','normalized','pos',[0 0 1 1]);
                case 'uicontrol'
                    L.hu = opt.in;
                    L.hp = get(opt.in,'parent');
                case 'axes'
                    % replace axes by an uicontrol
                    ha = opt.in;
                    L.hp = get(ha,'parent');
                    L.hu = uicontrol('parent',L.hp,'units',get(ha,'units'),'pos',get(ha,'pos'));
                    delete(ha)
                otherwise
                    error('bad handle')
            end
            L.hf = fn_parentfigure(L.hp);
            
            % context menu
            initlocalmenu(L)
            
            % list and event (bottom-up)
            set(L.hu,'style','listbox','min',0,'max',2, ...
                'callback',@(hu,evnt)event(L,'select'), ...
            	'keypressfcn',@(hu,evnt)keypress(L,evnt))
            if fn_switch(L.scrollwheel)
                L.scrollwheel = 'on'; % this will automaticall register scroll wheel
            end
            if isempty(get(L.hf,'WindowButtonMotionFcn'))
                % force update of current position when moving the mouse
                % around
                set(L.hf,'WindowButtonMotionFcn',@(u,e)donothing())
            end
            
            % watch filter
            connectlistener(F,L,'ChangedOperation',@(u,e)displayselection(L));
            connectlistener(F,L,'ChangedPoint',@(u,e)displaycross(L));
            
            % auto-delete
            set(L.hu,'deletefcn',@(u,evnt)delete(L))

            % update display (here, just sets the correct value)
            preformatvalues(L)
            displayselection(L)
            displaycross(L)
            displaylabel(L)
            
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
            set(L.hu,'UIContextMenu',m)
            
            uimenu(m,'label','new singleton selections','callback',@(u,e)event(L,'newuni'))
            uimenu(m,'label','new group selection [A]','callback',@(u,e)event(L,'newgroup'))
            uimenu(m,'label','add to selection','callback',@(u,e)event(L,'add'))

            uimenu(m,'label','define new group...','callback',@(u,e)event(L,'definegroup'),'separator','on')

            uimenu(m,'label','sort selections according to list order','callback',@(u,e)event(L,'sort'),'separator','on')
            uimenu(m,'label','reorder selections...','callback',@(u,e)event(L,'reorder'))
            
            uimenu(m,'label','remove highlighted group(s)','callback',@(u,e)event(L,'rmgroup'),'separator','on')
            uimenu(m,'label','remove all groups','callback',@(u,e)event(L,'rmgroupall'))
            uimenu(m,'label','remove highlighted individuals','callback',@(u,e)event(L,'rmuni'))
            uimenu(m,'label','remove all individuals','callback',@(u,e)event(L,'rmuniall'))
            uimenu(m,'label','remove highlighted selections','callback',@(u,e)event(L,'rm'))
            uimenu(m,'label','remove all selections','callback',@(u,e)event(L,'rmall'))

            L.menuitems.selmultin = uimenu(m,'separator','on','checked',fn_switch(L.selmultin), ...
                'label','temporary selection: individuals','callback',@(u,e)set(L,'selmultin',~L.selmultin));
            uimenu(m,'label','select all','callback',@(u,e)event(L,'selectall'))
            
            m1 = uimenu(m,'label','scroll wheel','separator','on');
            L.menuitems.scrollwheel = uimenu(m1,'label','activated', ...
                'checked',L.scrollwheel,'callback',@(u,e)set(L,'scrollwheel',fn_switch(L.scrollwheel,'toggle')));
            uimenu(m1,'label','make default in figure', ...
                'callback',@(u,e)set(L,'scrollwheel','default'));
        end
        function delete(L)
            delete@xplr.graphnode(L)
            if ~isvalid(L) && ~isprop(L,'hu'), return, end
            deleteValid(L.hu,L.hlabel)
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
            val = get(L.hu,'value');
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
                    set(L.hu,'value',1:L.F.szin)
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
                    set(L.hu,'value',val,'listboxtop',val(1)-2);
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
                        if iscell(val), newsel = val; else newsel = {val}; end
                        if L.doroi
                            for i=1:length(newsel)
                                newsel{i} = xplr.selectionnd('point1D',newsel{i});
                            end
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
                    if L.doroi
                        newsel = xplr.selectionnd.empty(1,0);
                    else
                        newsel = cell(1,0);
                    end
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
                updateSelection(L.F,'chg',idxrm,newsel,'SoftSelection',newissoft)
            elseif nnew>nrm
                updateSelection(L.F,'chg&new',{idxrm nsel+(1:nnew-nrm)},newsel,'SoftSelection',newissoft)
            elseif nrm>nnew
                updateSelection(L.F,'chg&rm',{idxrm(1:nnew) idxrm(nnew+1:nrm)},newsel,'SoftSelection',newissoft)
            end
        end
        function sel = buildCurrentSelection(L,domultin)
            val = get(L.hu,'value');
            if domultin
                if L.doroi
                    if L.F.headerin.ismeasure && all(diff(val)==1)
                        % selection is a segment rather than a mere
                        % list of points
                        sel = xplr.selectionnd('line1D',val([1 end]));
                    else
                        sel = xplr.selectionnd('point1D',val);
                    end
                else
                    sel = num2cell(val);
                end
            else
                if L.doroi
                    sel = xplr.selectionnd('point1D',val);
                else
                    sel = {val};
                end
            end
        end
    end
        
    % Get/Set - scroll wheel
    methods
        function set.scrollwheel(L,flag)
            switch flag
                case 'on'
                    fn_scrollwheelregister(L.hu,@(n)event(L,'scroll',n)) %#ok<MCSUP>
                    L.scrollwheel = 'on';
                case 'default'
                    fn_scrollwheelregister(L.hu,@(n)event(L,'scroll',n),'default') %#ok<MCSUP>
                    L.scrollwheel = 'on';
                case 'off'
                    fn_scrollwheelregister(L.hu,flag) %#ok<MCSUP>
                    L.scrollwheel = 'off';
                otherwise
                    error 'scrollwheel value must be ''off'', ''on'' or ''default'''
            end
            set(L.menuitems.scrollwheel,'checked',L.scrollwheel) %#ok<MCSUP>
        end
    end
    
    % Get/Set
    methods
        function set.selmultin(L,val)
            if val==L.selmultin, return, end
            L.selmultin = val;
            % update menu item
            set(L.menuitems.selmultin,'checked',fn_switch(val)) %#ok<MCSUP>
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
            top = get(L.hu,'listboxtop'); % the portion of the list that is shown moves when setting the string -> reset it to its current position
            set(L.hu,'string',str,'ListboxTop',top)
            if nsel==0
                % if no selection, highlight the point selection
                displaycross(L)
            else
                set(L.hu,'value',[selinds{softsel}])
            %             elseif isfield(L.F.shared,'list')
            %                 set(L.hu,'value',L.F.shared.list.cursel)
            end
        end
        function displaycross(L)
            set(L.hu,'value',L.F.index);
        end
        function displaylabel(L)
            fullfigure = isequal(fn_pixelsize(L.hu),fn_pixelsize(L.hp));
            if ~isempty(L.hlabel)
                % a control has already be created by the user for the
                % label, it only needs to be filled in
                set(L.hlabel,'style','text','string',L.F.headerout.label, ...
                    'horizontalalignment','center');
            elseif fullfigure
                % list occupies the full figure, put label in figure name
                % if parent is the figure
                if L.hp==L.hf, set(L.hf,'name',L.F.headerout.label), end
            else
                % list occupies only part of the figure, put label above
                % the list
                L.hlabel = uicontrol('parent',L.hp,'style','text','string',L.F.headerout.label, ...
                    'horizontalalignment','center');
                fn_controlpositions(L.hlabel,L.hu,[0 1 1 0],[0 2 0 15])
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
   