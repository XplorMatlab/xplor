classdef viewdisplay < hgsetget
    
    % Content
    properties (SetAccess='private')
        % data
        V           % parent 'view' object
        %L zoomslicer
        % graphics
        hp          % display panel
        ha          % main axes
        % display
        nodisplay = false   % data is too large, cancel display
        htransform  % containers for line/images that will be translated/scaled
        hdisplay    % handles of line/images
        %L1 labels      % xplrlight.displaylabel object
        graph       % xplrlight.displaygraph object
        %L navigation  % xplrlight.displaynavigation object
        %L2 clipping    % xplrlight.cliptool object
        %L3 % listeners that need being deleted upon object deletion
        %L3 listeners %= struct('slice',[],'zslice',[],'axsiz',[],'clip',[]); %L ,'zoom',[]
    end
    % Some "working memory"
    properties (Access='private')
        sliceChangeEvent
    end
    % Display properties
    properties (SetObservable=true, AbortSet=true)
        displaymode = 'time courses';  % 'time courses' or 'image' %L -> changed default to 'time courses' instead of 'images'
    end
    properties (SetAccess='private')
        org                                 % set with function setOrg
        activedim = struct('x',[],'y',[])   % change with function makeDimActive(D,d)
        colordim = [];                      % set with setColorDim
        clip = [0 1]                        % set with setClip, auto-clip with autoClip, other clip settings with sub-object cliptool
    end
    % Fast access (dependent)
    properties (Dependent, SetAccess='private')
        nd
        slice
        zslice
        %L zoomfilters
    end
    
    % Constructor, destructor
    methods
        function D = viewdisplay(V)
            % parent 'view' object and panel
            D.V = V;
            D.hp = V.panels.display;
            set(D.hp,'deletefcn',@(u,e)delete(D))
            
            %L % zoom slicer zooms into "slice" to yield "zslice"
            %L D.zoomslicer = xplrlight.zoomslicer(V.slicer.slice);

            % axes
            D.ha = axes('parent',D.hp);
            axis(D.ha,[-.5 .5 -.5 .5]) % center 0, available space 1
            set(D.ha,'box','on') %L4 ,'clim',[0 1])
            %L3 D.listeners.axsiz = fn_pixelsizelistener(D.ha,@(u,e)axisresize(D));
            %L3 c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position during all the following initializations
            
            %L4 % 'time courses'/'image' switch
            %L4 p = fn_propcontrol(D,'displaymode',{'popupmenu' 'time courses' 'image'},'parent',D.hp);
            %L4 fn_controlpositions(p.hu,D.hp,[1 1],[-90 -25 90 25])
            
            % positionning (needed by both labels and data display)
            D.graph = xplrlight.displaygraph(D);
            
            %L1 % automatic label positionning
            %L1 D.labels = xplrlight.displaylabels(D);
            
            %L % navigation (sliders, mouse actions)
            %L D.navigation = xplrlight.displaynavigation(D);
            
            %L2 % clipping tool
            %L2 D.clipping = xplrlight.cliptool(V.hf);
            %L2 D.listeners.clip = addlistener(D.clipping,'ChangedClip',@(u,e)clipchange(D,e));
            
            % set organization, connect sliders, display data and labels
            D.sliceChangeEvent = struct('flag','global');
            zslicechange(D)
            
            % listeners
            %L D.listeners.slice = addlistener(D.slice,'ChangedData',@(u,e)set(D,'sliceChangeEvent',e)); % mark that slice has changed, but treat it only later
            %L D.listeners.zoom = addlistener(D.zoomslicer,'ChangedZoom',@(u,e)zoomchange(D,e));
            %L3 D.listeners.zslice = addlistener(D.zslice,'ChangedData',@(u,e)zslicechange(D,e));
        end
        %L3 function delete(D)
        %L3     if ~isvalid(D) || ~isprop(D,'listeners'), return, end
        %L3     deleteValid(D.listeners,D.hp)
        %L3 end
    end
    
    % Dependent properties
    methods
        function n = get.nd(D)
            n = D.V.slicer.slice.nd;
        end
        function x = get.slice(D)
            x = D.V.slicer.slice;
        end
        function x = get.zslice(D)
            %L x = D.zoomslicer.slice;
            x = D.slice; % L
        end
        %L function zoomfilters = get.zoomfilters(D)
        %L     zoomfilters = [D.zoomslicer.filters.obj];
        %L end
    end
    
    % Change organization and active dim
    methods (Access='private')
        function out = checkActiveDim(D,doImmediateUpdate)
            % function anychg = checkActiveDim(D,doImmediateUpdate)
            %---
            % check that current active dims are ok, but also set active
            % dims if there aren't
            %L5 if nargin<2, doImmediateUpdate = true; end
            curorg = D.org;
            sz = D.slice.sz; okdim = (sz>1);
            if ~isempty(curorg.x), okdim(curorg.x(1)) = true; end
            if strcmp(D.displaymode,'image') && ~isempty(curorg.y), okdim(curorg.y(1)) = true; end
            
            % invalid dimensions
            % (x)
            dx = D.activedim.x;
            if ~isempty(dx) && ~(okdim(dx) && any(dx==[curorg.x curorg.yx]))
                dx = [];
            end
            % (y)
            dy = D.activedim.y;
            if ~isempty(dy) && ~(okdim(dy) && any(dy==[curorg.y curorg.xy]))
                dy = [];
            end
            
            % set new values
            if isempty([dx dy]) && ~isempty([curorg.xy curorg.yx])
                % (xy/yx)
                if ~isempty(curorg.xy)
                    dy = curorg.xy;
                else
                    dx = curorg.yx;
                end
            else
                % (x)
                if isempty(dx)
                    xorgok = curorg.x; xorgok(~okdim(xorgok)) = [];
                    if ~isempty(xorgok), dx = xorgok(end); end
                end
                % (y)
                if isempty(dy)
                    yorgok = curorg.y; yorgok(~okdim(yorgok)) = [];
                    if ~isempty(yorgok), dy = yorgok(end); end
                end
            end
            
            % update property
            newactd = struct('x',dx,'y',dy);
            anychg = ~isequal(newactd,D.activedim);
            if nargout>0, out = anychg; end
            if ~anychg, return, end
            D.activedim = newactd;
            
            %L5 % update display
            %L5 if doImmediateUpdate
            %L     D.navigation.connectFilter()
            %L1     D.labels.updateLabels('active')
            %L5     D.graph.setTicks()
            %L5 end
        end
    end
    methods
        function set.org(D,neworg)
            % check
            ok = isequal(fieldnames(neworg),{'x' 'y' 'ystatic' 'xy' 'yx'}') ...
                && ~(strcmp(D.displaymode,'image') && ~isempty(neworg.ystatic)) ...
                && length([neworg.xy neworg.yx])<=1;
            if ~ok, error 'organization structure is invalid', end
            % set
            D.org = neworg;
        end
        function setOrg(D,neworg,doImmediateDisplay)
            % function setOrg(D,neworg[,doImmediateDisplay])
            %---
            % if doImmediateDisplay is set to false, only labels are
            % updated; if it is set to true, update happens regardless of
            % whether neworg is actually new or not (this allows finishing
            % a previous incomplete update with doImmediateDisplay set to
            % false)
            %L3 c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            if nargin<3
                if isequal(neworg,D.org), return, end
                doImmediateDisplay = true;
            end
            D.org = neworg;
            % first update graph (new positionning will be needed for both
            % labels and data display)
            D.graph.computeSteps()
            % update labels
            if doImmediateDisplay, D.checkActiveDim(false), end
            D.labels.updateLabels()
            % update ticks and display
            if ~doImmediateDisplay, return, end
            D.checkColorDim(false)
            D.graph.setTicks()
            updateDisplay(D)
            % update slider connections
            connectFilter(D.navigation)
        end
        function makeDimActive(D,d,flag)
            %L3 c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            dotoggle = nargin>=3 && strcmp(flag,'toggle');
            % update active dim and connect slider
            if ismember(d,[D.org.x D.org.yx])
                if dotoggle && any(d==D.activedim.x)
                    D.activedim.x = []; 
                else
                    D.activedim.x = d;
                    if ismember(d,D.org.yx) || any(ismember(D.activedim.y,[D.org.xy D.org.yx]))
                        D.activedim.y = [];
                        D.navigation.connectFilter('y')
                    end
                end
                D.navigation.connectFilter('x')
            elseif ismember(d,[D.org.y D.org.xy])
                if dotoggle && any(d==D.activedim.y)
                    D.activedim.y = []; 
                else
                    D.activedim.y = d;
                    if ismember(d,D.org.xy) || any(ismember(D.activedim.x,[D.org.xy D.org.yx]))
                        D.activedim.x = [];
                        D.navigation.connectFilter('x')
                    end
                end
                D.navigation.connectFilter('y')
            end
            % update ticks and labels
            D.graph.setTicks()
            D.labels.updateLabels('active');
        end
    end
    
    %L7 % Color
    %L7 methods
    %L7     function setColorDim(D,d,doImmediateDisplay)
    %L7         if ~isnumeric(d) || numel(d)>1, error 'colordim must be empty or scalar', end
    %L7         if isequal(d,D.colordim), return, end
    %L7         if ~isempty(d) && ~isempty(D.org.x) && d==D.org.x(1), disp 'first x-dimension cannot be used for colors', return, end
    %L7         if nargin<3, doImmediateDisplay = true; end
    %L7         % set property
    %L7         D.colordim = d;
    %L7         % update display
    %L7         if doImmediateDisplay && strcmp(D.displaymode,'time courses')
    %L7             D.updateDisplay('color')
    %L7         end
    %L7     end
    %L7     function checkColorDim(D,doImmediateDisplay)
    %L7         cdim = D.colordim;
    %L7         if ~isempty(cdim) && ~isempty(D.org.x) && cdim==D.org.x(1)
    %L7             % cannot color according to the first x-dimension
    %L7             if nargin<2, doImmediateDisplay = false; end
    %L7             D.setColorDim([],doImmediateDisplay)
    %L7         end
    %L7     end
    %L7 end
    
    % Clipping
    methods
        function setClip(D,clip,doupdatedisplay)
            if ~isnumeric(clip) || length(clip)~=2 || diff(clip)<=0 || any(isnan(clip)|isinf(clip))
                disp 'clip value is not valid'
                return
            end
            if all(clip==D.clip), return, end
            if nargin<3, doupdatedisplay = true; end
            % set property
            D.clip = clip;
            % update display
            if doupdatedisplay, updateDisplay(D,'clip'), end
        end
        %L2 function autoClip(D,doupdatedisplay)
        %L2     if nargin<2, doupdatedisplay = true; end
        %L2     try
        %L2         val = fn_clip(D.zslice.data(:),D.clipping.autoclipmode,'getrange');
        %L2         setClip(D,val,doupdatedisplay)
        %L2     catch ME
        %L2         disp(ME)
        %L2     end
        %L2 end
        %L2 function clipchange(D,e)
        %L2     switch e.flag
        %L2         case 'clip'
        %L2             D.setClip(e.value)
        %L2         case 'automode'
        %L2             D.autoClip()
        %L2         case 'adjust'
        %L2             D.updateDisplay('clip')
        %L2         case 'span'
        %L2             if strcmp(D.clipping.span,'curview')
        %L2                 D.autoClip()
        %L2             end
        %L2         otherwise
        %L2             error('invalid ChangedClip flag ''%s''',e.flag)
        %L2     end
        %L2 end
    end
    
    % Update display
    methods (Access = 'private')
        function slicechange(D,e)
            % function slicechange(D,e)
            %---
            % function slicechange updates the 'org' property and slider
            % connections, but does not update display
            % it should be called after the zoomslicer has already
            % completed its update after slice change (hence the use of the
            % 'sliceChangeEvent' property to delay call to slicechange)
            
            if nargin<2, flag = 'global'; else flag = e.flag; end
            
            % Update organization
            switch flag
                case 'global'
                    % default organization: 1st and 4th dimension in x, 2nd
                    % and all other dimensions in y
                    neworg = struct('x',[],'y',[],'ystatic',[],'xy',[],'yx',[]);
                    nd = D.nd;
                    doimage = strcmp(D.displaymode,'image');
                    if nd==0, return, end
                    if nd==3 && doimage, neworg.xy = 3; end
                    if nd>=4, neworg.x = [1 4]; elseif nd>=1, neworg.x = 1; end
                    if nd>=3+doimage, neworg.y = [2 3 5:nd]; elseif nd>=2, neworg.y = 2; end
                    if ~isequal(neworg,D.org), D.org = neworg; end
                case 'insertdim'
                    error 'not implemented yet'
                case 'rmdim'
                    error 'not implemented yet'
                otherwise
                    % no change in organization, including in the 'chgdim'
                    % case!
            end
            
            % Update active dim and slider connections
            if fn_ismemberstr(flag,{'global' 'insertdim' 'rmdim' 'permdim'})
                D.checkActiveDim(false)
                %L D.navigation.connectFilter()
            elseif fn_ismemberstr(flag,{'chgdata' 'chg'})
                % slice size did not change
            else
                D.checkActiveDim(false)
                %L D.navigation.connectFilter()
            end
            
            %L7 % Update color dim
            %L7 D.checkColorDim(false)
        end
        function updateDisplay(D,flag,dim,ind)
            % function updateDisplay(D[,flag,dim,ind])
            if nargin<3, dim = []; end
            
            % Is data too large for being displayed?
            % TODO: this check should occur before setting the ticks
            %L sztest = D.zslice.sz;
            %L if ~isempty(D.org.x), sztest(D.org.x(1)) = 1; end
            %L if strcmp(D.displaymode,'image') && ~isempty(D.org.y), sztest(D.org.y(1)) = 1; end
            %L if strcmp(D.displaymode,'time courses')
            %L     ok = prod(sztest) <= xplrlight.parameters.get('display.NLineMax');
            %L else
            %L     ok = prod(sztest) <= xplrlight.parameters.get('display.NImageMax');
            %L end
            %L if ~ok
            %L     % too many grid elements: cancel display!
            %L     D.nodisplay = true;
            %L     deleteValid(D.htransform) % this will also delete children D.hdisplay
            %L     delete(findall(D.ha,'type','text','tag','xytick'))
            %L     set(D.ha,'xtick',[],'ytick',[])
            %L     if isempty(findall(D.ha,'type','text','tag','nodisplay'))
            %L         text(0,0,{'DATA IS TOO LARGE AND CANNOT BE DISPLAYED' 'BIN IT, OR USE FILTERS TO SLICE IT'}, ...
            %L             'parent',D.ha,'horizontalalignment','center','tag','nodisplay')
            %L     end
            %L     return
            %L elseif D.nodisplay
            %L     % display was canceled last time: we need a global update
            %L     D.nodisplay = false;
            %L     delete(findall(D.ha,'type','text','tag','nodisplay'))
            %L     flag = 'global';
            %L end                
            
            %L % Show watch
            %L c = fn_watch(D.V.hf); %#ok<NASGU>

            % To really run fast, avoid accessing object properties
            % repeatedly: access them once for all here
            dotimecourses = strcmp(D.displaymode,'time courses');
            %L2 clipadjust = D.clipping.adjust;
            %L2 if strcmp(clipadjust,'mean(line)')
            %L2     clipadjust = fn_switch(dotimecourses,'mean','none');
            %L2 end
            
            % What to do
            if nargin<2, flag = 'global'; end
            if ~strcmp(flag,'global')           %L5 
                error 'flag not handled yet'    %L5 
            end                                 %L5 
            %L5 if ~fn_ismemberstr(flag,{'clip' 'global' 'chgdata' 'chgdata&blocksize' 'new' 'remove' 'chg' 'perm' 'pos' 'color'})
            %L5     error 'flag not handled'
            %L5 end
            %L5 doreset = strcmp(flag,'global');
            %L5 donew = strcmp(flag,'new');
            %L5 doremove = strcmp(flag,'remove');
            %L5 dodispatch = ~fn_ismemberstr(flag,{'chg' 'chgdata' 'clip' 'color'});
            %L5 doselectdata = fn_ismemberstr(flag,{'chg' 'new'});
            %L5 doalldata = fn_ismemberstr(flag,{'global' 'chgdata' 'chgdata&blocksize' 'perm' 'clip' 'color'}); % color is set when updating ydata, but updating ydata is actually not necessary when only color changes...
            %L5 dochgx = strcmp(flag,'chgdata&blocksize');
            %L6 docolor = dotimecourses && ~fn_ismemberstr(flag,{'chgdata' 'chgdata&blocksize' 'clip'});
            
            % Reshape slice data adequately
            sz = D.zslice.sz;
            if ~isempty(D.org.x)
                xorg0 = D.org.x(1);
            else
                xorg0 = length(sz)+1;
                sz(end+1) = 1;
            end
            x = D.zslice.data;
            if dotimecourses
                nt = sz(xorg0);
                if nt==1
                    lineopt = {'linestyle','none','marker','.'};
                else
                    lineopt = {'linestyle','-','marker','none'};
                end
                szo = [prod(sz(1:xorg0-1)) prod(sz(xorg0+1:end))];
                x = reshape(x,[szo(1) nt szo(2)]);
                yorg0 = [];
            else
                %L4 if ~isempty(D.org.y)
                %L4     yorg0 = D.org.y(1);
                %L4 else
                %L4     yorg0 = length(sz)+1;
                %L4     sz(end+1) = 1;
                %L4 end
                %L4 dotranspose = (yorg0<xorg0);
                %L4 if dotranspose, [xorg0 yorg0] = deal(yorg0,xorg0); end
                %L4 nx = sz(xorg0);
                %L4 ny = sz(yorg0);
                %L4 szo = [prod(sz(1:xorg0-1)) prod(sz(xorg0+1:yorg0-1)) prod(sz(yorg0+1:end))];
                %L4 x = reshape(x,[szo(1) nx szo(2) ny szo(3)]);
            end
            
            %L5 % Check that current htransform and hdisplay are valid
            %L5 nd = D.zslice.nd;
            %L5 sz1 = sz; sz1([xorg0 yorg0]) = 1;
            %L5 sz1prev = sz1; 
            %L5 if ~isempty(dim), sz1prev(dim) = sz1prev(dim)+(doremove-donew)*length(ind); end
            %L5 if ~isequal(xplrlight.strictsize(D.htransform,nd),sz1prev) || ~all(ishandle(D.htransform(:)))...
            %L5         || ~isequal(xplrlight.strictsize(D.htransform,nd),sz1prev) || ~all(ishandle(D.htransform(:)))
            %L5     [doreset dodispatch doalldata] = deal(true);
            %L5     doselectdata = false;
            %L5 end
            
            %L6 % Prepare color
            %L6 if docolor
            %L6     cdim = D.colordim;
            %L6     colorhead = D.zslice.header(cdim);
            %L6     if isempty(D.colordim)
            %L6         if ~doreset, set(D.hdisplay(:),'color','k'), end
            %L6         docolor = false;
            %L6     else
            %L6         kcolor = strcmp({colorhead.sublabels.label},'ViewColor');
            %L6         if any(kcolor)
            %L6             cmap = cell2mat(colorhead.values(:,kcolor));
            %L6         else
            %L6             cmap = fn_colorset('plot12',1:colorhead.n);
            %L6         end
            %L6     end
            %L6 end
            
            % Prepare display and grid
            sz1 = sz; sz1([xorg0 yorg0]) = 1;
            %L5 if doreset          % reset display and grid elements
                deleteValid(D.htransform) % this will also delete children D.hdisplay
                [D.htransform D.hdisplay] = deal(gobjects([sz1 1]));
            %L5     [dodispatch doalldata] = deal(true);
            %L5 elseif donew      	% new grid elements
            %L5     subs = substruct('()',repmat({':'},1,D.zslice.nd));
            %L5     subs.subs{dim} = ind;
            %L5     D.htransform = subsasgn(D.htransform,subs,gobjects);
            %L5     D.hdisplay = subsasgn(D.hdisplay,subs,gobjects);
            %L5 elseif doremove     % remove grid elements
            %L5     subs = substruct('()',repmat({':'},1,D.zslice.nd));
            %L5     subs.subs{dim} = ind;
            %L5     deleteValid(subsref(D.htransform,subs)) % this also deletes the children hdisplay objects
            %L5     D.htransform = subsasgn(D.htransform,subs,[]);
            %L5     D.hdisplay = subsasgn(D.hdisplay,subs,[]);
            %L5 end
            
            % Prepare clipping
            clip0 = D.clip;
            clipextent = diff(clip0);
            %L2 if ~strcmp(clipadjust,'none'), clip0 = clip0-mean(clip0); end
            
            % Prepare several list of indices beforehand to avoid repeated
            % calls to functions such as ind2sub
            idxlist = 1:prod(sz1);
            %L5 if doselectdata && ~dodispatch
            %L5     % not all grid elements need to be visited ('chg' flag)
            %L5     idxlist = reshape(idxlist,sz1);
            %L5     subs = substruct('()',repmat({':'},1,D.zslice.nd));
            %L5     subs.subs{dim} = ind;
            %L5     idxlist = row(subsref(idxlist,subs));
            %L5 end
            ijklist = fn_indices(sz1,idxlist,'g2i');
            if dotimecourses
                [ibeflist iaftlist] = ind2sub(szo,idxlist);
            else
                [ibeflist imidlist iaftlist] = ind2sub(szo,idxlist);
            end
            
            % Prepare dispatch
            %L5 if dodispatch
                if dotimecourses
                    M = D.graph.gettransform(ijklist,[0 1]);
                else
                    M = D.graph.gettransform(ijklist);
                end
            %L5 end
            
            % Go! Loop on grid elements
            for u = 1:length(idxlist)
                idx = idxlist(u);
                %L5 ijk = ijklist(:,u);
                %L5 docurdata = doalldata || (doselectdata && any(ind==ijk(dim)));
                %L5 docreatecur = doreset || (donew && docurdata);
                % container
                %L5 if dodispatch
                %L5     if docreatecur
                        D.htransform(idx) = hgtransform('parent',D.ha,'matrix',M(:,:,idx),'HitTest','off');
                %L5     else
                %L5         set(D.htransform(idx),'matrix',M(:,:,idx))
                %L5     end
                %L5 end
                % line/image
                %L5 if docreatecur || docurdata
                    % get the data and adjust clipping if requested
                    if dotimecourses
                        xi = x(ibeflist(u),:,iaftlist(u));
                    else
                        xi = x(ibeflist(u),:,imidlist(u),:,iaftlist(u));
                    end
                    %L2 switch clipadjust
                    %L2     case 'none'
                    %L2         clipi = clip0;
                    %L2     case 'mean'
                    %L2         % adjustment by the mean
                    %L2         clipi = clip0 + nmean(xi(:));
                    %L2     otherwise
                    %L2         error('unknown clipping adjustment flag ''%s''',clipadjust)
                    %L2 end
                    clipi = clip0; %L2
                    % display it
                    if dotimecourses
                        xi = (xi-clipi(1))/clipextent;
                        %L5 if docreatecur
                            hl = line(1:nt,xi, ...
                                'parent',D.htransform(idx),'HitTest','off',lineopt{:});
                            D.hdisplay(idx) = hl;
                            %L6 if docolor, set(hl,'color',cmap(ijk(cdim),:)), end
                        %L5 else
                        %L5     hl = D.hdisplay(idx);
                        %L5     if dochgx
                        %L5         set(hl,'xdata',1:nt,lineopt{:})
                        %L5     end
                        %L5     set(hl,'ydata',xi)
                        %L5     %L6 if docolor, set(hl,'color',cmap(ijk(cdim),:)), end
                        %L5 end
                    else
                        %L4 xi = fn_clip(xi,clipi,[0 1]);
                        %L4 if dotranspose
                        %L4     im = permute(xi,[2 4 1 3 5]);
                        %L4 else
                        %L4     im = permute(xi,[4 2 1 3 5]);
                        %L4 end
                        %L4 if docreatecur
                        %L4     % y coordinates are negative to implement an up-down flip (D.graph.ystep(1) is negative,
                        %L4     % but the y-scale of the hgtransform cannot be set negative, so the absolute value is set and this hack on coordinates is used)
                        %L4     D.hdisplay(idx) = surface([.5 size(im,2)+.5],[-.5 -.5-size(im,1)],zeros(2), ...
                        %L4         'parent',D.htransform(idx), ...
                        %L4         'EdgeColor','none','FaceColor','texturemap','CDataMapping','scaled','CData',im, ...
                        %L4         'HitTest','off');
                        %L4 elseif dochgx
                        %L4     set(D.hdisplay(idx),'xdata',[.5 size(im,2)+.5],'ydata',[-.5 -.5-size(im,1)],'CData',im)
                        %L4 else
                        %L4     set(D.hdisplay(idx),'CData',im)
                        %L4 end
                    end
                %L5 end
            end
        end
    end
    methods
        function zslicechange(D,e)
            if nargin<2, flag = 'global'; else flag = e.flag; end
            %L3 c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            
            % Did slice change as well?
            if ~isempty(D.sliceChangeEvent)
                slicechange(D,D.sliceChangeEvent)
                D.sliceChangeEvent = [];
            end
            
            % Update graph (will be needed by both labels and data display)
            %L5 prevsz = D.graph.zslicesz;
            D.graph.computeSteps()
            
            % Update ticks and labels
            if ~(strcmp(flag,'chgdata') || (strcmp(flag,'chg') && ~any(e.dim==[D.activedim.x D.activedim.y])))
                D.graph.setTicks()
            end
            %L1 if fn_ismemberstr(flag,{'all' 'new' 'remove' 'chg&new' 'chg&rm' 'global' 'chgdim' 'insertdim' 'rmdim' 'permdim'})
            %L1     switch flag
            %L1         case 'global'
            %L1             D.labels.updateLabels('global')
            %L1         case {'chgdim' 'insertdim' 'rmdim' 'permdim'}
            %L1             D.labels.updateLabels(flag,e.dim)
            %L1         otherwise
            %L1             D.labels.updateLabels()
            %L1     end
            %L1 end
            
            % Update clipping
            chgclip = strcmp(flag,'global'); %L2 || strcmp(D.clipping.span,'curview');
            %L2 if chgclip, autoClip(D,false), end
            
            % Update clipping and display
            if fn_ismemberstr(flag,{'global' 'chgdim' 'insertdim' 'rmdir'})
                % Reset display
                updateDisplay(D,'global')
            else                                %L5 
                error('not implemented yet')    %L5 
            %L5 elseif strcmp(flag,'chgdata')
            %L5     % No change in size, all data need to be redisplayed
            %L5     updateDisplay(D,'chgdata')
            %L5 else
            %L5     % Smart display update
            %L5     dim = e.dim;
            %L5     if (~isempty(D.org.x) && D.org.x(1)==dim) ...
            %L5             || (strcmp(D.displaymode,'image') && ~isempty(D.org.y) && D.org.y(1)==dim)
            %L5         % changes are within elements (the grid arrangement
            %L5         % remains the same)
            %L5         if fn_ismemberstr(flag,{'perm' 'chg'}) ...
            %L5                 || (strcmp(flag,'all') && D.zslice.header(dim).n==prevsz(dim));
            %L5             flag = 'chgdata'; % no change in size
            %L5         else
            %L5             flag = 'chgdata&blocksize';
            %L5         end
            %L5         updateDisplay(D,flag)
            %L5     elseif ~chgclip
            %L5         % the grid arrangement changes
            %L5         switch flag
            %L5             case 'chg'  % check this case first because it is this one that occurs when going fast through a list
            %L5                 updateDisplay(D,'chg',e.dim,e.ind)
            %L5             case {'new' 'remove' 'perm'}
            %L5                 updateDisplay(D,flag,dim,e.ind)
            %L5             case 'all'
            %L5                 ncur = size(D.htransform,dim);
            %L5                 n = D.zslice.sz(dim);
            %L5                 if n==ncur
            %L5                     updateDisplay(D,'chgdata')
            %L5                 elseif n>ncur
            %L5                     updateDisplay(D,'new',dim,ncur+1:n)
            %L5                     updateDisplay(D,'chg',dim,1:ncur)
            %L5                 else
            %L5                     updateDisplay(D,'remove',dim,n+1:ncur)
            %L5                     updateDisplay(D,'chg',dim,1:n)
            %L5                 end
            %L5             case 'chg&new'
            %L5                 updateDisplay(D,'new',e.dim,e.ind{2})
            %L5                 updateDisplay(D,'chg',e.dim,e.ind{1})
            %L5             case 'chg&rm'
            %L5                 updateDisplay(D,'remove',e.dim,e.ind{2})
            %L5                 updateDisplay(D,'chg',e.dim,e.ind{1})
            %L5             otherwise
            %L5                 error('flag ''%s'' is not handled','flag')
            %L5         end
            %L5     elseif chgclip
            %L5         % all grid elements need to be updated
            %L5         ncur = size(D.htransform,dim);
            %L5         n = D.zslice.sz(dim);
            %L5         if n==ncur
            %L5             updateDisplay(D,'chgdata')
            %L5         elseif n>ncur
            %L5             updateDisplay(D,'chgdata')
            %L5             updateDisplay(D,'new',dim,ncur+1:n)
            %L5         else
            %L5             updateDisplay(D,'remove',dim,n+1:ncur)
            %L5             updateDisplay(D,'chgdata')
            %L5         end
            %L5     end
            end
        end
        %L function zoomchange(D,e)
        %L     % update graph positions: if data has changed in size,
        %L     % positioning will be updated upon notification of data change;
        %L     % however if data has not changed in size, positioning needs to
        %L     % be updated here
        %L     if ~e.chgnout
        %L         c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
        %L         D.graph.computeSteps()
        %L         D.graph.setTicks()
        %L         D.labels.updateLabels()
        %L         updateDisplay(D,'pos')
        %L     end
        %L end
        %L3 function axisresize(D)
        %L3     c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
        %L3     D.graph.computeSteps()
        %L3     D.graph.setTicks()
        %L3     updateDisplay(D,'pos')
        %L3 end
        %L4 function set.displaymode(D,mode)
        %L4     %L3 c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            % set property
        %L4     D.displaymode = mode;
        %L4     % for 'image' mode, check that org is valid, and modify it if
        %L4     % necessary
        %L4     if strcmp(mode,'image') && ~isempty(D.org.ystatic)
        %L4         % it is not possible to superimpose images -> change org
        %L4         neworg = D.org;
        %L4         neworg.y = [D.org.y D.org.ystatic];
        %L4         neworg.ystatic = [];
        %L4         D.setOrg(neworg) % this automatically updates display among other things
        %L4     else
        %L4         % update display
        %L4         D.graph.computeSteps()
        %L4         D.graph.setTicks()
        %L4         D.labels.updateLabels()
        %L4         updateDisplay(D)
        %L4     end
        %L4 end
        %L function resetDisplay(D)
        %L     % reset axis
        %L     cla(D.ha)
        %L     % re-display everything
        %L     zslicechange(D)
        %L end
    end
end