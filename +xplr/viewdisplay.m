classdef viewdisplay < hgsetget
    
    % Content
    properties (SetAccess='private')
        % data
        V           % parent 'view' object
        zoomslicer
        % graphics
        hp          % display panel
        ha          % main axes
        % display
        nodisplay = false   % data is too large, cancel display
        htransform  % containers for line/images that will be translated/scaled
        hdisplay    % handles of line/images
        hlegend     % handle of legend axes
        labels      % xplr.displaylabel object
        graph       % xplr.displaygraph object
        navigation  % xplr.displaynavigation object
        clipping    % xplr.cliptool object
        colormap    % xplr.colormap object
        % listeners that need being deleted upon object deletion
        listeners %= struct('slice',[],'zslice',[],'zoom',[],'axsiz',[],'clip',[]);
    end
    % Some "working memory"
    properties (Access='private')
        sliceChangeEvent
    end
    % Display properties
    properties (SetObservable=true, AbortSet=true)
        displaymode = 'image';  % 'time courses' or 'image'
        showcolorlegend = false; % false by default because of bug in Matlab's legend function, which inactivates several listeners
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
        zoomfilters
    end
    
    % Constructor, destructor
    methods
        function D = viewdisplay(V)
            % parent 'view' object and panel
            D.V = V;
            D.hp = V.panels.display;
            set(D.hp,'deletefcn',@(u,e)delete(V))
            
            % zoom slicer zooms into "slice" to yield "zslice"
            D.zoomslicer = xplr.zoomslicer(V.slicer.slice,D);

            % axes
            D.ha = axes('parent',D.hp);
            axis(D.ha,[-.5 .5 -.5 .5]) % center 0, available space 1
            set(D.ha,'box','on','clim',[0 1])
            try set(D.ha, 'XTickLabelRotation',45), end % recent Matlab versions only
            try set(D.ha,'TickLabelInterpreter','none'), end % recent Matlab versions only
            D.listeners.axsiz = fn_pixelsizelistener(D.ha,@(u,e)axisresize(D));
            c = disableListener(D.listeners.axsiz); % prevent display update following automatic change of axis position during all the following initializations

            % 'time courses'/'image' switch
            p = fn_propcontrol(D,'displaymode',{'popupmenu' 'time courses' 'image'},'parent',D.hp);
            fn_controlpositions(p.hu,D.hp,[1 1],[-90 -25 90 25])
            
            % positionning (needed by both labels and data display)
            D.graph = xplr.displaygraph(D);
            
            % automatic label positionning
            D.labels = xplr.displaylabels(D);
            
            % navigation (sliders, mouse actions)
            D.navigation = xplr.displaynavigation(D);
            
            % clipping tool
            D.clipping = xplr.cliptool(V.hf); % creates a menu
            D.listeners.clip = addlistener(D.clipping,'ChangedClip',@(u,e)clipchange(D,e));
            
            % colormap tool
            D.colormap = xplr.colormaptool(D); % creates a menu
            D.listeners.colormap = addlistener(D.colormap,'ChangedColorMap',@(u,e)colormap(V.hf,D.colormap.cmap)); %#ok<CPROP>
            
            % set organization, connect sliders, display data and labels
            D.sliceChangeEvent = struct('flag','global');
            zslicechange(D)
            
            % listeners
            D.listeners.slice = addlistener(D.slice,'ChangedData',@(u,e)set(D,'sliceChangeEvent',e)); % mark that slice has changed, but treat it only later
            D.listeners.zoom = addlistener(D.zoomslicer,'ChangedZoom',@(u,e)zoomchange(D,e));
            D.listeners.zslice = addlistener(D.zslice,'ChangedData',@(u,e)zslicechange(D,e));
            
            % problem: c won't be deleted automatically (and axsiz listener
            % might not be re-enabled) because the workspace continue to
            % exist, because of all the anonymous functions that were
            % defined
            delete(c)
        end
        function delete(D)
            if ~isprop(D,'listeners'), return, end
            deleteValid(D.listeners,D.hp)
        end
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
            x = D.zoomslicer.slice;
        end
        function zoomfilters = get.zoomfilters(D)
            zoomfilters = [D.zoomslicer.filters.obj];
        end
    end
    
    % Some usefull simple methods
    methods
        function s = getSize(D, unit, dim)
            % function s = getSize(D, unit [, dim])
            s = fn_objectsize(D.ha, unit);
            if nargin >= 3
                if isnumeric(dim)
                    s = s(dim);
                elseif strcmp(dim, 'x')
                    s = s(1);
                elseif strcmp(dim, 'y')
                    s = s(2);
                else
                    error argument
                end
            end
        end
    end
    
    % Change organization and active dim
    methods (Access='private')
        function out = checkActiveDim(D,doImmediateUpdate)
            % function anychg = checkActiveDim(D,doImmediateUpdate)
            %---
            % check that current active dims are ok, but also set active
            % dims if there aren't
            if nargin<2, doImmediateUpdate = true; end
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
            
            % update display
            if doImmediateUpdate
                D.navigation.connectFilter()
                D.labels.updateLabels('active')
                D.graph.setTicks()
            end
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
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            if nargin<3
                if isequal(neworg,D.org), return, end
                doImmediateDisplay = true;
            end
            D.org = neworg;
            % is zslice too large for being displayed
            D.checkzslicesize()
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
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
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
    
    % Color
    methods
        function setColorDim(D,d,doImmediateDisplay)
            if ~isnumeric(d) || numel(d)>1, error 'colordim must be empty or scalar', end
            if isequal(d,D.colordim), return, end
            if ~isempty(d) && ~isempty(D.org.x) && d==D.org.x(1), disp 'first x-dimension cannot be used for colors', return, end
            if nargin<3, doImmediateDisplay = true; end
            % set property
            D.colordim = d;
            % update color legend?
            if D.showcolorlegend
                displayColorLegend(D)
            end
            % update display
            if doImmediateDisplay && strcmp(D.displaymode,'time courses')
                D.updateDisplay('color')
            end
        end
        function checkColorDim(D,doImmediateDisplay)
            cdim = D.colordim;
            if ~isempty(cdim) && ~isempty(D.org.x) && cdim==D.org.x(1)
                % cannot color according to the first x-dimension
                if nargin<2, doImmediateDisplay = false; end
                D.setColorDim([],doImmediateDisplay)
            end
        end
        function set.showcolorlegend(D,val)
            D.showcolorlegend = val;
            displayColorLegend(D)
        end
        function displayColorLegend(D)
            % delete current legend
            delete(D.hlegend)
            D.hlegend = [];
            % do really display a legend
            d = D.colordim;
            if isempty(d) || ~D.showcolorlegend || strcmp(D.displaymode,'image'), return, end
            % get line handles
            s = substruct('()',num2cell(ones(1,D.nd)));
            s.subs{d} = ':';
            hl = subsref(D.hdisplay,s);
            % display legend
            names = D.slice.header(d).getItemNames;
            D.hlegend = fn_colorlegend(row(hl),names,'SouthWest','frame');
        end
    end
    
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
        function autoClip(D,doupdatedisplay)
            if nargin<2, doupdatedisplay = true; end
            try
                val = fn_clip(D.zslice.data(:),D.clipping.autoclipmode,'getrange');
                setClip(D,val,doupdatedisplay)
            catch ME
                disp(ME)
            end
        end
        function clipchange(D,e)
            switch e.flag
                case 'clip'
                    D.setClip(e.value)
                case 'automode'
                    D.autoClip()
                case 'adjust'
                    D.updateDisplay('clip')
                case 'span'
                    if strcmp(D.clipping.span,'curview')
                        D.autoClip()
                    end
                otherwise
                    error('invalid ChangedClip flag ''%s''',e.flag)
            end
        end
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
                D.navigation.connectFilter()
            elseif fn_ismemberstr(flag,{'chgdata' 'chg'})
                % slice size did not change
            else
                D.checkActiveDim(false)
                D.navigation.connectFilter()
            end
            
            % Update color dim
            D.checkColorDim(false)
        end
        function updateDisplay(D,flag,dim,ind)
            % function updateDisplay(D[,flag,dim,ind])
            if nargin<3, dim = []; end
            
            % Is data too large for being displayed?
            if D.nodisplay
                % too many grid elements: cancel display!
                deleteValid(D.htransform) % this will also delete children D.hdisplay
                delete(findall(D.ha,'type','text','tag','xytick'))
                set(D.ha,'xtick',[],'ytick',[])
                if isempty(findall(D.ha,'type','text','tag','nodisplay'))
                    text(0,0,{'DATA IS TOO LARGE AND CANNOT BE DISPLAYED' 'BIN IT, OR USE FILTERS TO SLICE IT'}, ...
                        'parent',D.ha,'horizontalalignment','center','tag','nodisplay')
                end
                return
            end
            hnodisplay = findall(D.ha,'type','text','tag','nodisplay');
            if ~isempty(hnodisplay)
                % display was canceled last time: we need a global update
                delete(hnodisplay)
                flag = 'global';
            end                
            
            % Show watch
            c = fn_watch(D.V.hf); %#ok<NASGU>

            % To really run fast, avoid accessing object properties
            % repeatedly: access them once for all here
            dotimecourses = strcmp(D.displaymode,'time courses');
            clipadjust = D.clipping.adjust;
            if strcmp(clipadjust,'mean(line)')
                clipadjust = fn_switch(dotimecourses,'mean','none');
            end
            
            % What to do
            if nargin<2, flag = 'global'; end
            if ~fn_ismemberstr(flag,{'clip' 'global' 'chgdata' 'chgdata&blocksize' 'new' 'remove' 'chg' 'perm' 'pos' 'color'})
                error 'flag not handled'
            end
            doreset = strcmp(flag,'global');
            donew = strcmp(flag,'new');
            doremove = strcmp(flag,'remove');
            dodispatch = ~fn_ismemberstr(flag,{'chg' 'chgdata' 'clip' 'color'});
            doselectdata = fn_ismemberstr(flag,{'chg' 'new'});
            doalldata = fn_ismemberstr(flag,{'global' 'chgdata' 'chgdata&blocksize' 'perm' 'clip' 'color'}); % color is set when updating ydata, but updating ydata is actually not necessary when only color changes...
            dochgx = strcmp(flag,'chgdata&blocksize');
            docolor = dotimecourses && ~fn_ismemberstr(flag,{'chgdata' 'chgdata&blocksize' 'clip'});
            
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
                if ~isempty(D.org.y)
                    yorg0 = D.org.y(1);
                else
                    yorg0 = length(sz)+1;
                    sz(end+1) = 1;
                end
                dotranspose = (yorg0<xorg0);
                if dotranspose, [xorg0 yorg0] = deal(yorg0,xorg0); end
                nx = sz(xorg0);
                ny = sz(yorg0);
                szo = [prod(sz(1:xorg0-1)) prod(sz(xorg0+1:yorg0-1)) prod(sz(yorg0+1:end))];
                x = reshape(x,[szo(1) nx szo(2) ny szo(3)]);
            end
            
            % Check that current htransform and hdisplay are valid
            nd = D.zslice.nd;
            sz1 = sz; sz1([xorg0 yorg0]) = 1;
            sz1prev = sz1; 
            if ~isempty(dim), sz1prev(dim) = sz1prev(dim)+(doremove-donew)*length(ind); end
            if ~isequal(xplr.strictsize(D.htransform,nd),sz1prev) || ~all(ishandle(D.htransform(:)))...
                    || ~isequal(xplr.strictsize(D.htransform,nd),sz1prev) || ~all(ishandle(D.htransform(:)))
                [doreset dodispatch doalldata] = deal(true);
                doselectdata = false;
            end
            
            % Prepare color
            if docolor
                cdim = D.colordim;
                colorhead = D.zslice.header(cdim);
                if isempty(D.colordim)
                    if ~doreset, set(D.hdisplay(:),'color','k'), end
                    docolor = false;
                else
                    kcolor = strcmp({colorhead.sublabels.label},'ViewColor');
                    if any(kcolor)
                        cmap = cell2mat(colorhead.values(:,kcolor));
                    else
                        cmap = fn_colorset('plot12',1:colorhead.n);
                    end
                end
            end
            
            % Prepare display and grid
            sz1 = sz; sz1([xorg0 yorg0]) = 1;
            if doreset          % reset display and grid elements
                deleteValid(D.htransform) % this will also delete children D.hdisplay
                [D.htransform D.hdisplay] = deal(gobjects([sz1 1]));
                [dodispatch doalldata] = deal(true);
            elseif donew      	% new grid elements
                subs = substruct('()',repmat({':'},1,D.zslice.nd));
                subs.subs{dim} = ind;
                D.htransform = subsasgn(D.htransform,subs,gobjects);
                D.hdisplay = subsasgn(D.hdisplay,subs,gobjects);
            elseif doremove     % remove grid elements
                subs = substruct('()',repmat({':'},1,D.zslice.nd));
                subs.subs{dim} = ind;
                deleteValid(subsref(D.htransform,subs)) % this also deletes the children hdisplay objects
                D.htransform = subsasgn(D.htransform,subs,[]);
                D.hdisplay = subsasgn(D.hdisplay,subs,[]);
            end
            
            % Prepare clipping
            clip0 = D.clip;
            clipextent = diff(clip0);
            if ~strcmp(clipadjust,'none'), clip0 = clip0-mean(clip0); end
            
            % Prepare several list of indices beforehand to avoid repeated
            % calls to functions such as ind2sub
            idxlist = 1:prod(sz1);
            if doselectdata && ~dodispatch
                % not all grid elements need to be visited ('chg' flag)
                idxlist = reshape(idxlist,sz1);
                subs = substruct('()',repmat({':'},1,D.zslice.nd));
                subs.subs{dim} = ind;
                idxlist = row(subsref(idxlist,subs));
            end
            ijklist = fn_indices(sz1,idxlist,'g2i');
            if dotimecourses
                [ibeflist iaftlist] = ind2sub(szo,idxlist);
            else
                [ibeflist imidlist iaftlist] = ind2sub(szo,idxlist);
            end
            
            % Prepare dispatch
            if dodispatch
                if dotimecourses
                    M = D.graph.gettransform(ijklist,[0 1]);
                else
                    M = D.graph.gettransform(ijklist);
                end
            end
            
            % Go! Loop on grid elements
            for u = 1:length(idxlist)
                idx = idxlist(u);
                ijk = ijklist(:,u);
                docurdata = doalldata || (doselectdata && any(ind==ijk(dim)));
                docreatecur = doreset || (donew && docurdata);
                % container
                if dodispatch
                    if docreatecur
                        D.htransform(idx) = hgtransform('parent',D.ha,'matrix',M(:,:,idx),'HitTest','off');
                    else
                        set(D.htransform(idx),'matrix',M(:,:,idx))
                    end
                end
                % line/image
                if docreatecur || docurdata
                    % get the data and adjust clipping if requested
                    if dotimecourses
                        xi = x(ibeflist(u),:,iaftlist(u));
                    else
                        xi = x(ibeflist(u),:,imidlist(u),:,iaftlist(u));
                    end
                    switch clipadjust
                        case 'none'
                            clipi = clip0;
                        case 'mean'
                            % adjustment by the mean
                            clipi = clip0 + nmean(xi(:));
                        otherwise
                            error('unknown clipping adjustment flag ''%s''',clipadjust)
                    end
                    xi = double(xi);
                    clipi = double(clipi);
                    % display it
                    if dotimecourses
                        xi = (xi-clipi(1))/clipextent;
                        if docreatecur
                            hl = line(1:nt,xi, ...
                                'parent',D.htransform(idx),'HitTest','off',lineopt{:});
                            D.hdisplay(idx) = hl;
                            if docolor, set(hl,'color',cmap(ijk(cdim),:)), end
                        else
                            hl = D.hdisplay(idx);
                            if dochgx
                                set(hl,'xdata',1:nt,lineopt{:})
                            end
                            set(hl,'ydata',xi)
                            if docolor, set(hl,'color',cmap(ijk(cdim),:)), end
                        end
                    else
                        xi = fn_clip(xi,clipi,[0 1]);
                        if dotranspose
                            im = permute(xi,[2 4 1 3 5]);
                        else
                            im = permute(xi,[4 2 1 3 5]);
                        end
                        if docreatecur
                            % y coordinates are negative to orient the
                            % image downward (see also comment inside of
                            % displaygaph.gettransform method, where the
                            % y-scale of the hgtransform cannot be set
                            % negative)
                            D.hdisplay(idx) = surface([.5 size(im,2)+.5],[-.5 -.5-size(im,1)],zeros(2), ...
                                'parent',D.htransform(idx), ...
                                'EdgeColor','none','FaceColor','texturemap','CDataMapping','scaled','CData',im, ...
                                'HitTest','off');
                        elseif dochgx
                            set(D.hdisplay(idx),'xdata',[.5 size(im,2)+.5],'ydata',[-.5 -.5-size(im,1)],'CData',im)
                        else
                            set(D.hdisplay(idx),'CData',im)
                        end
                    end
                end
            end
        end
        function checkzslicesize(D)
            sztest = D.zslice.sz;
            if ~isempty(D.org.x), sztest(D.org.x(1)) = 1; end
            if strcmp(D.displaymode,'image') && ~isempty(D.org.y), sztest(D.org.y(1)) = 1; end
            if strcmp(D.displaymode,'time courses')
                ok = (prod(sztest) <= xplr.parameters.get('display.NLineMax')) && ...
                    (prod(D.zslice.sz) <= xplr.parameters.get('display.NLinePointMax'));
            else
                ok = (prod(sztest) <= xplr.parameters.get('display.NImageMax')) && ...
                    (prod(D.zslice.sz) <= xplr.parameters.get('display.NImagePixelMax'));
            end
            D.nodisplay = ~ok;
        end
    end
    methods
        function zslicechange(D,e)
            if nargin<2, flag = 'global'; else flag = e.flag; end
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            
            % Did slice change as well?
            if ~isempty(D.sliceChangeEvent)
                slicechange(D,D.sliceChangeEvent)
                D.sliceChangeEvent = [];
            end
            
            % Is zslice too large for being displayed
            D.checkzslicesize()
            
            % Update graph (will be needed by both labels and data display)
            prevsz = D.graph.zslicesz;
            D.graph.computeSteps()
            
            % Update ticks and labels
            if ~(strcmp(flag,'chgdata') || (strcmp(flag,'chg') && ~any(e.dim==[D.activedim.x D.activedim.y])))
                D.graph.setTicks()
            end
            if fn_ismemberstr(flag,{'all' 'new' 'remove' 'chg&new' 'chg&rm' 'global' 'chgdim' 'insertdim' 'rmdim' 'permdim'})
                switch flag
                    case 'global'
                        D.labels.updateLabels('global')
                    case {'chgdim' 'insertdim' 'rmdim' 'permdim'}
                        D.labels.updateLabels(flag,e.dim)
                    otherwise
                        D.labels.updateLabels()
                end
            end
            
            % Update clipping
            chgclip = strcmp(flag,'global') || strcmp(D.clipping.span,'curview');
            if chgclip, autoClip(D,false), end
            
            % Update clipping and display
            if fn_ismemberstr(flag,{'global' 'chgdim' 'insertdim' 'rmdir'})
                % Reset display
                updateDisplay(D,'global')
            elseif strcmp(flag,'chgdata')
                % No change in size, all data need to be redisplayed
                updateDisplay(D,'chgdata')
            else
                % Smart display update
                dim = e.dim;
                if (~isempty(D.org.x) && D.org.x(1)==dim) ...
                        || (strcmp(D.displaymode,'image') && ~isempty(D.org.y) && D.org.y(1)==dim)
                    % changes are within elements (the grid arrangement
                    % remains the same)
                    if fn_ismemberstr(flag,{'perm' 'chg'}) ...
                            || (strcmp(flag,'all') && D.zslice.header(dim).n==prevsz(dim));
                        flag = 'chgdata'; % no change in size
                    else
                        flag = 'chgdata&blocksize';
                    end
                    updateDisplay(D,flag)
                elseif ~chgclip
                    % the grid arrangement changes
                    switch flag
                        case 'chg'  % check this case first because it is this one that occurs when going fast through a list
                            updateDisplay(D,'chg',e.dim,e.ind)
                        case {'new' 'remove' 'perm'}
                            updateDisplay(D,flag,dim,e.ind)
                        case 'all'
                            ncur = size(D.htransform,dim);
                            n = D.zslice.sz(dim);
                            if n==ncur
                                updateDisplay(D,'chgdata')
                            elseif n>ncur
                                updateDisplay(D,'new',dim,ncur+1:n)
                                updateDisplay(D,'chg',dim,1:ncur)
                            else
                                updateDisplay(D,'remove',dim,n+1:ncur)
                                updateDisplay(D,'chg',dim,1:n)
                            end
                        case 'chg&new'
                            updateDisplay(D,'new',e.dim,e.ind{2})
                            updateDisplay(D,'chg',e.dim,e.ind{1})
                        case 'chg&rm'
                            updateDisplay(D,'remove',e.dim,e.ind{2})
                            updateDisplay(D,'chg',e.dim,e.ind{1})
                        otherwise
                            error('flag ''%s'' is not handled','flag')
                    end
                elseif chgclip
                    % all grid elements need to be updated
                    ncur = size(D.htransform,dim);
                    n = D.zslice.sz(dim);
                    if n==ncur
                        updateDisplay(D,'chgdata')
                    elseif n>ncur
                        updateDisplay(D,'chgdata')
                        updateDisplay(D,'new',dim,ncur+1:n)
                    else
                        updateDisplay(D,'remove',dim,n+1:ncur)
                        updateDisplay(D,'chgdata')
                    end
                end
            end
            
            % Update legend
            if fn_ismemberstr(flag,{'global' 'chgdim' 'insertdim' 'rmdir'})
                D.colordim = [];
            elseif ~strcmp(flag,'chgdata') && isequal(e.dim,D.colordim)
                displayColorLegend(D)
            end
        end
        function zoomchange(D,e)
            % update graph positions: if data has changed in size,
            % positioning will be updated upon notification of data change;
            % however if data has not changed in size, positioning needs to
            % be updated here
            if ~e.chgnout
                c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
                D.checkzslicesize % is zslice too large for being displayed
                D.graph.computeSteps()
                D.graph.setTicks()
                D.labels.updateLabels()
                updateDisplay(D,'pos')
            end
        end
        function axisresize(D)
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            D.graph.computeSteps()
            D.graph.setTicks()
            updateDisplay(D,'pos')
        end
        function set.displaymode(D,mode)
            c = disableListener(D.listeners.axsiz); %#ok<MCSUP,NASGU> % prevent display update following automatic change of axis position
            % set property
            D.displaymode = mode;
            % for 'image' mode, check that org is valid, and modify it if
            % necessary
            if strcmp(mode,'image') && ~isempty(D.org.ystatic)
                % it is not possible to superimpose images -> change org
                neworg = D.org;
                neworg.y = [D.org.y D.org.ystatic];
                neworg.ystatic = [];
                D.setOrg(neworg) % this automatically updates display among other things
            else
                % update display
                D.checkzslicesize() % is zslice too large for being displayed
                D.graph.computeSteps() %#ok<MCSUP>
                D.graph.setTicks() %#ok<MCSUP>
                D.labels.updateLabels() %#ok<MCSUP>
                updateDisplay(D)
            end
            % show/hide color legend
            displayColorLegend(D)
        end
        function resetDisplay(D)
            % reset axis
            cla(D.ha)
            % re-display everything
            zslicechange(D)
        end
    end
end