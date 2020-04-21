classdef viewdisplay < xplr.graphnode
    % viewdisplay
    
    % Content
    properties (SetAccess='private')
        % data
        V           % parent 'view' object
        zoomslicer
        previousheaders = xplr.header.empty(1,0);
        % graphics
        hp          % display panel
        ha          % main axes
        menu        % general display menu
        % display
        nodisplay = false   % data is too large, cancel display
        htransform  % containers for line/images that will be translated/scaled
        hdisplay    % handles of line/images
        gridclip    % clipping for each grid element
        hlegend     % handle of legend axes
        labels      % xplr.displaylabel object
        graph       % xplr.displaygraph object
        navigation  % xplr.displaynavigation object
        clipping    % xplr.cliptool object
        colormap    % xplr.colormap object
    end
    
    % Some "working memory"
    properties (Access='private')
        sliceChangeEvent
        listeners = struct;
    end
    
    % Display properties
    properties (SetObservable=true, AbortSet=true)
        displaymode = 'image';  % 'time courses' or 'image'
        showcolorlegend = false; % false by default because of bug in Matlab's legend function, which inactivates several listeners
        linealpha = 1; % lines have a small degree of transparency!
    end    
    properties (SetAccess='private')
        layoutID                              % layout, i.e. which data dimension appear on which location; set with function setLayoutID
        layoutIDall                           % layout, including singleton dimensions; set with function setLayoutID
        activedimID = struct('x',[],'y',[])   % dimensions on which zooming mouse actions and sliders apply; change with function makeDimActive(D,d)
        colordimID = [];                      % set with setColorDim
        clip = [0 1]                          % set with setClip, auto-clip with autoClip, other clip settings with sub-object cliptool
    end
    
    % Shortcuts (dependent)
    properties (Dependent, SetAccess='private')
        nd
        slice
        zslice
        zoomfilters
        activedim
        colordim
        layout
        internal_dim
        internal_dimID
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
            D.addListener(D.ha,D.listeners.axsiz);
            c = disableListener(D.listeners.axsiz); % prevent display update following automatic change of axis position during all the following initializations
            
            % 'time courses'/'image' switch
            p = fn_propcontrol(D,'displaymode',{'popupmenu' 'time courses' 'image'},'parent',D.hp);
            fn_controlpositions(p.hu,D.hp,[1 1],[-90 -25 90 25])
            
            % positionning (needed by both labels and data display)
            D.graph = D.addComponent(xplr.displaygraph(D));
            
            % automatic label positionning
            D.labels = D.addComponent(xplr.displaylabels(D));
            
            % general display menu
            D.menu = uimenu(D.V.hf,'label','Display','callback',@(u,e)display_menu(D));
            
            % clipping tool
            D.clipping = D.addComponent(xplr.cliptool(V.hf)); % creates a menu
            D.addListener(D.clipping,'ChangedClip',@(u,e)clipchange(D,e));
            
            % colormap tool
            D.colormap = D.addComponent(xplr.colormaptool(D)); % creates a menu
            D.addListener(D.colormap,'ChangedColorMap',@(u,e)D.updateDisplay('clip'));
            
            % navigation (sliders, mouse actions)
            D.navigation = D.addComponent(xplr.displaynavigation(D)); % creates a menu
            
            % set organization, connect sliders, display data and labels
            D.sliceChangeEvent = struct('flag','global');
            zslicechange(D)
            
            % listeners
            D.addListener(D.slice,'ChangedData',@(u,e)set(D,'sliceChangeEvent',e)); % mark that slice has changed, but treat it only later
            D.addListener(D.zoomslicer,'ChangedZoom',@(u,e)zoomchange(D,e));
            D.addListener(D.zslice,'ChangedData',@(u,e)zslicechange(D,e));
            
            % problem: c won't be deleted automatically (and axsiz listener
            % might not be re-enabled) because the workspace continue to
            % exist, because of all the anonymous functions that were
            % defined
            delete(c)
        end
        function delete(D)
            delete@xplr.graphnode(D)
            delete(D.zoomslicer)
            delete(D.navigation)
            delete(D.labels)
            delete(D.graph)
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
        function activedim = get.activedim(D)
            activedim = struct( ...
                'x', D.slice.dimensionNumber(D.activedimID.x), ...
                'y', D.slice.dimensionNumber(D.activedimID.y));                
        end
        function colordim = get.colordim(D)
            colordim = D.slice.dimensionNumber(D.colordimID);
            if isequal(colordim,0), colordim = []; end
        end
        function layout = get.layout(D)
            layout = D.layoutID.dimensionNumber();
        end
        function dim = get.internal_dim(D)
            org = D.layout;
            if strcmp(D.displaymode,'image') && ~isempty(org.y)
                if isempty(org.x)
                    dim = org.y(1);
                else
                    dim = [org.x(1) org.y(1)];
                end
            else
                if isempty(org.x)
                    dim = [];
                else
                    dim = org.x(1);
                end
            end
        end
        function dim = get.internal_dimID(D)
            org = D.layoutID;
            if strcmp(D.displaymode,'image') && ~isempty(org.y)
                if isempty(org.x)
                    dim = org.y(1);
                else
                    dim = [org.x(1) org.y(1)];
                end
            else
                if isempty(org.x)
                    dim = [];
                else
                    dim = org.x(1);
                end
            end
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
    
    % General methods
    methods (Access='private')
        function display_menu(D)
            m = D.menu;
            delete(get(m,'children'))
            
            % time courses display
            if strcmp(D.displaymode,'time courses')
                fn_propcontrol(D,'linealpha', ...
                    {'menuval', {1 .7 .4 .1}, {'none' 'mild' 'medium' 'strong' 'manual'}}, ...
                    {'parent',m,'label','Lines transparency'});
                dosep = true;
            else
                dosep = false;
            end
            
            % cross
            fn_propcontrol(D.navigation,'showcross','menu', ...
                {'parent',m,'label','Show cross','separator',onoff(dosep)});
            if D.navigation.showcross
                fn_propcontrol(D.navigation,'crosscolor', ...
                    {'menu', {'k' 'b' 'r' [1 1 1]*.6 'w'}, {'black' 'blue' 'red' 'gray' 'white' 'other'}}, ...
                    {'parent',m,'label','Cross color'});
                fn_propcontrol(D.navigation,'crossalpha', ...
                    {'menu', {1 .4 .05}, {'none' 'medium' 'barely visible' 'manual'}}, ...
                    {'parent',m,'label','Cross transparency'});
            end
            
            % separation marks
            org = D.layoutID;
            if length(org.x)>1 || length(org.y)>strcmp(D.displaymode,'image') || ~isempty([org.xy org.yx])
                fn_propcontrol(D.graph,'showseparation','menu', ...
                    {'parent',m,'label','Show separations between lines/images','separator','on'});
                if D.graph.showseparation
                    fn_propcontrol(D.graph,'separationcolor', ...
                        {'menu', {'k' [.8 .8 1] [1 .8 .8] [.8 .8 .8]}, {'black' 'light blue' 'light red' 'light gray' 'other'}}, ...
                        {'parent',m,'label','Separations color'});
                end
            end            
            
            % reset display
            uimenu(m,'label','Reset display','separator','on', ...
                'callback',@(u,e)D.resetDisplay())
        end
    end
    
    % Change filters, zoomfilters binning, organization and active dim
    methods (Access='private')
        function out = checkActiveDim(D,doImmediateUpdate,doAuto)
            % function anychg = checkActiveDim(D,doImmediateUpdate[,doAuto])
            %---
            % check that current active dims are ok, but also set active
            % dims if there aren't
            if nargin<2, doImmediateUpdate = true; end
            if nargin<3, doAuto = false; end
            
            orgID = D.layoutID;
            sz = D.slice.sz; 
            
            % try to keep same active dims if they still exist in the new
            % data
            dx = D.activedimID.x;
            if ~ismember(dx, [D.slice.header.dimID]), dx = []; end
            dy = D.activedimID.y;
            if ~ismember(dy, [D.slice.header.dimID]), dy = []; end
            
            % check position of active dims
            if isempty([dx dy]) 
                % no valid active dim: set some if doAuto is set to true
                if ~doAuto
                    % do not set active dims
                %                 elseif ~isempty(orgID.xy)
                %                     dy = orgID.xy;
                %                 elseif ~isempty(orgID.yx)
                %                     dx = orgID.yx;
                else
                    % (x)
                    if isempty(dx) && ~isempty(orgID.x)
                        dx = orgID.x(end);
                    end
                    % (y)
                    if isempty(dy) && ~isempty(orgID.y)
                        dy = orgID.y(end);
                    end
                end
            elseif ismember(dx, orgID.x)
                % ok, we can keep dx as the active x dimension
                if ismember(dy, orgID.y)
                    % we can also keep dy as the active y dimension
                    [dx, dy] = deal(dx, dy);
                else
                    [dx, dy] = deal(dx, []);
                end
            elseif ismember(dx, orgID.y)
                % dimension dx moved to y location
                if ismember(dy, orgID.x)
                    % and dimension dy moved to x location: switch the 2
                    % active dims!
                    [dx, dy] = deal(dy, dx);
                elseif ismember(dy, orgID.y)
                    % keep dy as the active y dimension, remove dx from
                    % active x
                    [dx, dy] = deal([], dy);
                else
                    % replace active y dimension with dx
                    [dx, dy] = deal([], dx);
                end
            else
                if ismember(dy, orgID.x)
                    [dx, dy] = deal(dy, []);
                elseif ismember(dy, orgID.y)
                    [dx, dy] = deal([], dy);
                elseif ismember(dx, orgID.xy)
                    [dx, dy] = deal([], dx);
                elseif ismember(dx, orgID.yx)
                    [dx, dy] = deal(dx, []);
                elseif ismember(dy, orgID.xy)
                    [dx, dy] = deal([], dy);
                elseif ismember(dy, orgID.yx)
                    [dx, dy] = deal(dy, []);
                else
                    % should not happen
                    [dx, dy] = deal([], []);
                end
            end
                               
            % update property
            newvalue = struct('x',dx,'y',dy);
            anychg = ~isequal(newvalue, D.activedimID);
            if nargout>0, out = anychg; end
            if ~anychg, return, end
            D.activedimID = newvalue;
            
            % update display
            if doImmediateUpdate
                D.navigation.connectZoomFilter()
                D.labels.updateLabels('active')
                D.graph.s()
                D.graph.setValueTicks()
            end
        end
    end    
    methods
        function dimensionContextMenu(D,m,dim)
            % function dimensionContextMenu(D,m,dim)
            %---
            % This function populates the context menu that appears when
            % right-clicking on a label
            
            [dim, dimID] = D.slice.dimensionNumberAndID(dim);
            head = D.slice.header(dim);
            delete(get(m,'children'))
            
            % Line properties
            % (color)
            docolor = strcmp(D.displaymode,'time courses');
            if docolor
                uimenu(m,'label',['Color according to ' head.label],'checked',onoff(isequal(D.colordimID,dimID)), ...
                    'callback',@(u,e)D.setColorDim(fn_switch(isequal(D.colordimID,dimID),[],dimID)))
                uimenu(m,'label','Display color legend', ...
                    'enable',onoff(isequal(D.colordimID,dimID)),'checked',onoff(D.showcolorlegend), ...
                    'callback',@(u,e)set(D,'showcolorlegend',~D.showcolorlegend))
            end

            % Binning
            m1 = uimenu(m,'label','Binning','Separator',onoff(docolor));
            binvalues = {1 2 3 4 5 10 20 'set'};
            bindisplays = {'none' '2' '3' '4' '5' '10' '20' 'other...'};
            curbin = D.zoomfilters(dim).bin;
            for i=1:length(binvalues)
                bin = binvalues{i};
                uimenu(m1,'label',bindisplays{i},'checked',onoff(isequal(curbin,bin)), ...
                    'callback',@(u,e)setbin(D,dim,bin));
            end

            % select ZoomFilter key (check the created menu item
            % that corresponds to the current key)
            m2 = uimenu(m,'label','zoom filter','Separator',onoff(docolor));
            availablekeys = xplr.bank.availableFilterKeys('zoomfilter');
            newkey = max(availablekeys)+1;
            keyvalues = [0 availablekeys newkey];
            fn_num2str(availablekeys, 'shared zoom %i', 'cell');
            keydisplays = [ ...
                'private zoom' ...
                fn_num2str(availablekeys, 'shared zoom %i', 'cell') ...
                num2str(newkey,'shared zoom %i (new key)')
                ];
            curkey = D.zoomfilters(dim).linkkey;
            for i=1:length(keyvalues)
                keyvalue = keyvalues(i);      
                uimenu(m2,'label',keydisplays{i},'checked',onoff(isequal(curkey,keyvalue)), ...
                    'callback',@(u,e)D.zoomslicer.changeKey(dim,keyvalue));
            end

            % select crossSelector key 
            curfilt = D.navigation.pointfilters{dim};
            if ~isempty(curfilt)
                m2 = uimenu(m,'label','cross selector key','Separator',onoff(docolor));

                availablekeys = xplr.bank.availableFilterKeys('point');
                newkey = max(availablekeys)+1;
                keyvalues = [0 availablekeys newkey];
                fn_num2str(availablekeys, 'cross selector key %i', 'cell');
                keydisplays = [ ...
                    'private cross selector' ...
                    fn_num2str(availablekeys, 'cross selector key %i', 'cell') ...
                    num2str(newkey,'cross selector key %i (new key)')
                    ];
                    curkey = curfilt.linkkey;
                uimenu(m2, 'label', 'show point selector','callback',@(u,e)xplr.bank.showList(curfilt));
                for i=1:length(keyvalues)
                    keyvalue = keyvalues(i);
                    uimenu(m2,'label',keydisplays{i},'checked',onoff(isequal(curkey,keyvalue)), ...
                        'callback',@(u,e)connectPointFilter(D.navigation,dim,keyvalue));
                end
            end
        end
        function setbin(D,d,bin)
            if strcmp(bin,'set')
                bin = fn_input('Binning',D.zoomfilters(d).bin,'stepper 1 1 Inf 1');
                if isempty(bin), return, end
            end
            D.zoomfilters(d).setBin(bin)
        end
        function setLayoutID(D,newlayoutID,doImmediateDisplay)
            % function setLayoutID(D,newlayoutID[,doImmediateDisplay])
            %---
            % if doImmediateDisplay is set to false, only labels are
            % updated; if it is set to true, update happens regardless of
            % whether newlayout is actually new or not (this allows finishing
            % a previous incomplete update with doImmediateDisplay set to
            % false)
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            if nargin<3
                if isequal(newlayoutID,D.layoutIDmem), return, end
                doImmediateDisplay = true;
            end
            D.layoutIDall = newlayoutID;
            D.layoutID = newlayoutID.currentlayout(); % keep only dimensions actually displayed
            % is zslice too large for being displayed
            D.checkzslicesize()
            % first update graph (new positionning will be needed for both
            % labels and data display)
            D.graph.computeSteps()
            % update labels
            if doImmediateDisplay, D.checkActiveDim(false), end
            D.labels.updateLabels()
            drawnow
            % check whether color dim and active dim remain valid
            D.checkColorDim(false)
            D.checkActiveDim(false,true)
            % update ticks and display
            if ~doImmediateDisplay, return, end
            D.graph.setTicks()
            updateDisplay(D) % will call setValueTicks if necessary
            % update slider connections
            connectZoomFilter(D.navigation)                      
            % reposition cross
            D.navigation.repositionCross()
            % update selection display
            D.navigation.displayselection()
        end
        function set_dim_location(D,dimID,location,doImmediateDisplay)
            % function set_dim_location(D,dimID,location,doImmediateDisplay)
            %---
            % set new location of specific dimensions; locations of other
            % dimensions will automatically be adjusted
            % more details in the help of xplr.displaylayout.set_dim_location
            %
            % See also xplr.displaylayout.set_dim_location
            
            % new layout
            newlayoutID = D.layoutIDall.set_dim_location(dimID,location);
            
            % update display
            if nargin<4, doImmediateDisplay = true; end
            D.setLayoutID(newlayoutID,doImmediateDisplay)
        end
        function set.activedimID(D,value)
            D.activedimID = value;
        end
        function makeDimActive(D,dimID,flag)
            dimID = D.slice.dimensionID(dimID);
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            dotoggle = nargin>=3 && strcmp(flag,'toggle');
            % update active dim and connect slider
            if ismember(dimID,[D.layoutID.x D.layoutID.yx])
                if dotoggle && any(dimID==D.activedimID.x)
                    D.activedimID.x = [];
                else
                    D.activedimID.x = dimID;
                    if ismember(dimID,D.layoutID.yx) || any(ismember(D.activedimID.y,[D.layoutID.xy D.layoutID.yx]))
                        D.activedimID.y = [];
                        D.navigation.connectZoomFilter('y')
                    end
                end
                D.navigation.connectZoomFilter('x')
            elseif ismember(dimID,[D.layoutID.y D.layoutID.xy])
                if dotoggle && any(dimID==D.activedimID.y)
                    D.activedimID.y = [];
                else
                    D.activedimID.y = dimID;
                    if ismember(dimID,D.layoutID.xy) || any(ismember(D.activedimID.x,[D.layoutID.xy D.layoutID.yx]))
                        D.activedimID.x = [];
                        D.navigation.connectZoomFilter('x')
                    end
                end
                D.navigation.connectZoomFilter('y')
            end
            % update ticks and labels
            D.graph.setTicks()
            D.graph.setValueTicks()
            D.labels.updateLabels('active');
        end
    end
    
    % Color
    methods
        function setColorDim(D,dim,doImmediateDisplay)
            if ~isnumeric(dim) || numel(dim)>1, error 'colordim must be empty or scalar', end
            dimID = D.slice.dimensionID(dim);
            if isequal(dimID,D.colordimID), return, end
            if ~isempty(dimID) && ~isempty(D.layoutID.x) && dimID==D.layoutID.x(1), disp 'first x-dimension cannot be used for colors', return, end
            if nargin<3, doImmediateDisplay = true; end
            % set property
            D.colordimID = dimID;
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
            cdimID = D.colordimID;
            if ~isempty(cdimID) && ~isempty(D.layoutID.x) && cdimID==D.layoutID.x(1)
                % cannot color according to the first x-dimension
                if nargin<2, doImmediateDisplay = false; end
                D.setColorDim([],doImmediateDisplay)
            end
        end
        function set.linealpha(D,linealpha)
            % check
            if ~isscalar(linealpha) || linealpha<=0 || linealpha>1, error 'incorrect alpha value for lines', end
            % set value
            D.linealpha = linealpha;
            % update display
            if strcmp(D.displaymode,'time courses')
                D.updateDisplay('color')
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
                xplr.debuginfo('stop','clip value is not valid')
                return
            end
            if all(clip==D.clip), return, end
            if nargin<3, doupdatedisplay = true; end
            % set property
            D.clip = double(clip);
            % update display
            if doupdatedisplay, updateDisplay(D,'clip'), end
        end
        function autoClip(D,doupdatedisplay)
            if nargin<2, doupdatedisplay = true; end
            try
                val = fn_clip(D.zslice.data(:),D.clipping.autoclipmode,'getrange');
                if isinf(val(1)), val(1) = -1e6; end
                if isinf(val(2)), val(2) = 1e6; end
                if ~any(isnan(val)), setClip(D,val,doupdatedisplay), end
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
    methods (Static)
        function ok = testDisplayable(sz,displaymode,layout)
            szgrid = sz;
            if ~isempty(layout.x), szgrid(layout.x(1)) = 1; end
            if strcmp(displaymode,'image') && ~isempty(layout.y), szgrid(layout.y(1)) = 1; end
            if strcmp(displaymode,'time courses')
                ok = (prod(szgrid) <= xplr.parameters.get('display.NLineMax')) && ...
                    (prod(sz) <= xplr.parameters.get('display.NLinePointMax'));
            else
                ok = (prod(szgrid) <= xplr.parameters.get('display.NImageMax')) && ...
                    (prod(sz) <= xplr.parameters.get('display.NImagePixelMax'));
            end
        end
    end
    methods (Access = 'private')
        function slicechange(D,e)
            % function slicechange(D,e)
            %---
            % function slicechange updates the 'layout' property and slider
            % connections, but does not update display
            % it should be called after the zoomslicer has already
            % completed its update after slice change (hence the use of the
            % 'sliceChangeEvent' property to delay call to slicechange)
            
            if nargin<2, flag = 'global'; else, flag = e.flag; end
            xplr.debuginfo('viewdisplay','slicechange %s', flag)
            
            % first time?
            if isempty(D.layoutIDall)
                % some heuristics to choose initial layout
                D.displaymode = fn_switch(sum(D.slice.sz>1) == 1, 'time courses', 'image');
                D.layoutIDall = xplr.displaylayout(D);
                D.layoutID = D.layoutIDall;
            else
                % keep locations of dimensions already present in
                % D.layoutIDall, use some heuristic to choose
                % locations of new dimensions
                [D.layoutIDall, D.layoutID] = D.layoutIDall.updateLayout();
            end
            
            % Update active dim and slider connections
            if fn_ismemberstr(flag,{'global'})
                D.checkActiveDim(false,true)
                D.navigation.connectZoomFilter()
            elseif fn_ismemberstr(flag,{'chgdata' 'chg'})
                % slice size did not change
            else
                D.checkActiveDim(false)
                D.navigation.connectZoomFilter()
            end
            
            % Update color dim
            D.checkColorDim(false)
            
            % Assign point filters to each updated dimension
            switch flag
                case 'chgdata'
                    % nothing to do: only the data has changed
                case 'global'
                    D.navigation.connectPointFilter()
                case {'all' 'chgdim' 'new' 'remove' 'chg' 'chg&new' 'chg&rm' 'perm'}
                    dim = D.slice.dimensionNumber(e.dim);
                    D.navigation.connectPointFilter(dim)
                otherwise
                    error('flag ''%s'' not handled', flag)
            end
            
            % Check whether current dimension for selections display is
            % still valid, i.e. whether the connected filter still fits the
            % dimension in the new slice (if it is still valid, note that
            % selection display update will occur in D.zslicechange)
            D.navigation.checkselectionfilter()
            
            % Se previous headers to current headers
            D.previousheaders = D.slice.header;
        end
        function updateDisplay(D,flag,dim,ind)
            % function updateDisplay(D[,flag,dim,ind])
            if nargin<3, dim = []; end
            
            % Is data too large for being displayed?
            if D.nodisplay
                % too many grid elements: cancel display!
                deleteValid(D.htransform) % this will also delete children D.hdisplay
                D.gridclip = [];
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
            
            % What to do
            if nargin<2, flag = 'global'; end
            if ~fn_ismemberstr(flag,{'clip' 'global' 'chgdata' 'chgdata&blocksize' 'new' 'remove' 'chg' 'perm' 'pos' 'color'})
                error 'flag not handled'
            end
            doreset = strcmp(flag,'global');
            donew = strcmp(flag,'new');
            doremove = strcmp(flag,'remove');
            doposition = ~fn_ismemberstr(flag,{'chg' 'chgdata' 'clip' 'color'});
            dodataall = fn_ismemberstr(flag,{'clip' 'global' 'chgdata' 'chgdata&blocksize' 'perm' 'color'}); % color is set when updating ydata, but updating ydata is actually not necessary when only color changes...
            dodataselect = fn_ismemberstr(flag,{'new' 'chg'});
            dodata = dodataall || dodataselect;
            dochgx = strcmp(flag,'chgdata&blocksize');
            docolor = dotimecourses && ~fn_ismemberstr(flag,{'chgdata' 'chgdata&blocksize' 'clip'});
            
            % Reshape slice data adequately: after reshape, dimensions in x
            % will alternate between external and internal dimensions
            sz = D.zslice.sz;
            org = D.layout; % convert from dim ID to dim numbers
            if ~isempty(org.x)
                xlayout0 = D.layout.x(1);
            else
                % no dimension displayed in x; mimick an additional
                % dimension displayed in x
                xlayout0 = length(sz)+1;
                sz(end+1) = 1;
            end
            % build a 3-element vector of 'internal' dimensions; some of
            % them might be fake dimensions, the idea is to have in any
            % case 3 internal dimensions to avoid handling multiple
            % separate cases
            if dotimecourses
                internaldim = [xlayout0 length(sz)+(1:2)];
                sz(end+(1:2)) = 1;
            else
                if ~isempty(org.y)
                    ylayout0 = org.y(1);
                else
                    ylayout0 = length(sz)+1;
                    sz(end+1) = 1;
                end
                if ~isempty(org.mergeddata)
                    channeldim = org.mergeddata;
                else
                    channeldim = length(sz)+1;
                    sz(end+1) = 1;
                end
                internaldim = [xlayout0 ylayout0 channeldim];
            end
            % prepare the dimensions permutation after extracting each
            % sub-signal/sub-image
            [internaldim_sort, perm_xi2slice] = sort(internaldim);
            perm_slice2xi = [0 0 0 1 3 5 7];
            perm_slice2xi(perm_xi2slice) = [2 4 6];
            if ~dotimecourses
                % for image, XPLOR dimension order is x-y, but Matlab is
                % y-x, so we permute the two before displaying images
                perm_slice2xi(1:2) = perm_slice2xi([2 1]);
            end
            % size of reshape
            szr = zeros(1,7);
            for i = 1:3
                if i==1
                    szr(2*i-1) = prod(sz(1:internaldim_sort(1)-1));
                else
                    szr(2*i-1) = prod(sz(internaldim_sort(i-1)+1:internaldim_sort(i)-1));
                end
                szr(2*i) = sz(internaldim_sort(i));
            end
            szr(7) = prod(sz(internaldim_sort(3)+1:end));
            % reshape
            x = reshape(D.zslice.data, szr);
            % size of remaining external dimensions
            szo = szr(1:2:end);

            % Check that current htransform and hdisplay are valid
            nd = D.zslice.nd;
            sz1 = sz; sz1(internaldim) = 1;
            sz1prev = sz1;
            if ~isempty(dim), sz1prev(dim) = sz1prev(dim)+(doremove-donew)*length(ind); end
            if ~isequal(strictsize(D.htransform,nd),sz1prev) || ~all(ishandle(D.htransform(:)))...
                    || ~isequal(strictsize(D.htransform,nd),sz1prev) || ~all(ishandle(D.htransform(:)))
                [doreset doposition dodataall] = deal(true);
                dodataselect = false;
            end
            
            % Prepare color
            if docolor
                cdim = D.colordim;
                colorhead = D.zslice.header(cdim);
                if isempty(D.colordim)
                    if ~doreset, set(D.hdisplay(:),'color',[0 0 0 D.linealpha]), end
                    docolor = false;
                else
                    kcolor = strcmp({colorhead.sublabels.label},'ViewColor');
                    if any(kcolor)
                        cmap = cell2mat(colorhead.values(:,kcolor));
                    else
                        cmap = fn_colorset('plot12',1:colorhead.n);
                    end
                    if size(cmap,2) == 3 && D.linealpha < 1
                        cmap(:,4) = D.linealpha;
                    end
                end
            end
            
            % Prepare display and grid
            sz1 = sz; sz1(internaldim) = 1;
            if doreset          % reset display and grid elements
                deleteValid(D.htransform) % this will also delete children D.hdisplay
                [D.htransform, D.hdisplay] = deal(gobjects([sz1 1]));
                D.gridclip = zeros([2 sz1 1]);
                [doposition, dodataall] = deal(true);
            elseif donew      	% new grid elements
                subs = substruct('()',repmat({':'},1,D.zslice.nd));
                subs.subs{dim} = ind;
                D.htransform = subsasgn(D.htransform,subs,gobjects);
                D.hdisplay = subsasgn(D.hdisplay,subs,gobjects);
                subs.subs = [{':'} subs.subs];
                D.gridclip = subsasgn(D.gridclip,subs,0);
            elseif doremove     % remove grid elements
                subs = substruct('()',repmat({':'},1,D.zslice.nd));
                subs.subs{dim} = ind;
                deleteValid(subsref(D.htransform,subs)) % this also deletes the children hdisplay objects
                D.htransform = subsasgn(D.htransform,subs,[]);
                D.hdisplay = subsasgn(D.hdisplay,subs,[]);
                subs.subs = [{':'} subs.subs];
                D.gridclip = subsasgn(D.gridclip,subs,[]);
            end
            
            % Prepare clipping
            clip0 = D.clip;
            clipextent = diff(clip0);
            
            % Prepare several list of indices beforehand to avoid repeated
            % calls to functions such as ind2sub
            idxlist = 1:prod(sz1);
            if dodataselect && ~doposition
                % not all grid elements need to be visited ('chg' flag)
                idxlist = reshape(idxlist,sz1);
                subs = substruct('()',repmat({':'},1,D.zslice.nd));
                subs.subs{dim} = ind;
                idxlist = row(subsref(idxlist,subs));
            end
            ijklist = fn_indices(sz1,idxlist,'g2i');
            [i1, i2, i3, i4] = ind2sub(szo,idxlist);
            
            % Prepare dispatch
            if doposition
                M = D.graph.gettransform(ijklist);
            end
            
            % Go! Loop on grid elements
            for u = 1:length(idxlist)
                idx = idxlist(u);
                ijk = ijklist(:,u);
                dodatacur = dodataall || (dodataselect && any(ind==ijk(dim)));
                docreatecur = doreset || (donew && dodatacur);
                % container
                if doposition
                    if docreatecur
                        D.htransform(idx) = hgtransform('parent',D.ha,'matrix',M(:,:,idx),'HitTest','off');
                    else
                        set(D.htransform(idx),'matrix',M(:,:,idx))
                    end
                end
                % line/image
                if docreatecur || dodatacur
                    % get the data and adjust clipping if requested
                    xi = x(i1(u),:,i2(u),:,i3(u),:,i4(u));
                    xi = permute(xi, perm_slice2xi);
                    switch clipadjust
                        case 'none'
                            clipi = clip0;
                        case 'mean'
                            % adjustment by the mean
                            clipi = nmean(xi(:)) + [-.5 +.5] * clipextent;
                        otherwise
                            error('unknown clipping adjustment flag ''%s''',clipadjust)
                    end
                    xi = fn_float(xi);
                    clipi = clipi;
                    % store clipping values
                    D.gridclip(:,idx) = clipi;                    
                    % display it
                    if dotimecourses
                        xi = (xi-clipi(1))/clipextent;
                        nt = sz(internaldim(1));
                        if nt==1
                            lineopt = {'linestyle','none','marker','.'};
                        else
                            lineopt = {'linestyle','-','marker','none'};
                        end
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
                        nanvalue = .95;
                        % size in color dimension must be 1, 3 or 4;
                        % correct if it is not the case
                        alpha = [];
                        nc = sz(internaldim(3));
                        if nc == 2
                            % add a third blue channel set to zero
                            xi(:,:,3) = 0;
                        elseif nc > 4
                            xi = xi(:,:,1:3);
                        end
                        if size(xi,3) == 1
                            im = fn_clip(xi,clipi,D.colormap.cmap,nanvalue);
                        else
                            im = fn_clip(xi,clipi,[0 1],nanvalue);
                        end
                        if nc == 4
                            [im, alpha] = deal(im(:,:,1:3), im(:,:,4));
                        end
                        if docreatecur
                            % y coordinates are negative to orient the
                            % image downward (see also comment inside of
                            % displaygaph.gettransform method, where the
                            % y-scale of the hgtransform cannot be set
                            % negative)
                            D.hdisplay(idx) = surface([.5 size(im,2)+.5],[-.5 -.5-size(im,1)],zeros(2), ...
                                'parent',D.htransform(idx), ...
                                'EdgeColor','none','FaceColor','texturemap', ...
                                'CDataMapping','scaled', 'CData',im, ...
                                'HitTest','off');
                        elseif dochgx
                            set(D.hdisplay(idx),'xdata',[.5 size(im,2)+.5],'ydata',[-.5 -.5-size(im,1)],'CData',im)
                        else
                            set(D.hdisplay(idx),'CData',im)
                        end
                        if ~isempty(alpha)
                            set(D.hdisplay(idx),'FaceAlpha','texturemap','AlphaData',alpha)
                        end
                    end
                end
            end
            
            % update value y-ticks
            if dodata
                D.graph.setValueTicks()
            end
            
            % make sur containers are below labels, selections, etc.
            % (fater to have only one call to uistack rather than after
            % creating each element)
            if doreset || donew
                uistack(D.htransform(:), 'bottom')
            end

        end
        function checkzslicesize(D)
            D.nodisplay = ~D.testDisplayable(D.zslice.sz,D.displaymode,D.layout);
        end
    end
    methods
        function resetDisplay(D)
            % reset axis
            cla(D.ha)
            % re-display everything
%             D.sliceChangeEvent = struct('flag','global');
            D.navigation.displayselection('reset')
            zslicechange(D) % this will automatically re-create the cross, but not the selection displays
        end
    end
    
    % Callbacks (i.e. handle specific changes in slice, properties, etc.)
    methods
        function zslicechange(D,e)
            if nargin<2, flag = 'global'; else flag = e.flag; end
            xplr.debuginfo('viewdisplay','zslicechange %s',flag)
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            
            % Did slice change as well?
            if ~isempty(D.sliceChangeEvent)
                slicechange(D,D.sliceChangeEvent)
                D.sliceChangeEvent = [];
            end
            
            % Dimension(s) where change occured
            if nargin>=2
                [chgdim, chgdimID] = D.zslice.dimensionNumberAndID(e.dim);
            end
            
            % Is zslice too large for being displayed
            D.checkzslicesize()
            
            % Update graph (will be needed by both labels and data display)
            prevsz = D.graph.zslicesz;
            D.graph.computeSteps()
            
            % Update labels and ticks (do the labels first because ticks
            % update can change the size of the axes, and therefore trigger
            % labels re-positionning, which can cause error if the number
            % of labels has decreased)
            if fn_ismemberstr(flag,{'all' 'new' 'remove' 'chg&new' 'chg&rm' 'global' 'chgdim'})
                switch flag
                    case 'global'
                        D.labels.updateLabels('global')
                    case 'chgdim'
                        D.labels.updateLabels(flag,chgdim)
                    otherwise
                        D.labels.updateLabels()
                end
            end
            if ~(strcmp(flag,'chgdata') || (strcmp(flag,'chg') && ~any(chgdim==[D.activedim.x D.activedim.y])))
                D.graph.setTicks()
            end
            
            % Update clipping
            chgclip = strcmp(flag,'global') || strcmp(D.clipping.span,'curview');
            if chgclip, autoClip(D,false), end
            
            % Update display
            if fn_ismemberstr(flag,{'global' 'chgdim'})
                % Reset display
                updateDisplay(D,'global')
            elseif strcmp(flag,'chgdata')
                % No change in size, all data need to be redisplayed
                updateDisplay(D,'chgdata')
            else
                % Smart display update
                if (~isempty(D.layoutID.x) && D.layoutID.x(1)==chgdimID) ...
                        || (strcmp(D.displaymode,'image') && ~isempty(D.layoutID.y) && D.layoutID.y(1)==chgdimID)
                    % changes are within elements (the grid arrangement
                    % remains the same)
                    if fn_ismemberstr(flag,{'perm' 'chg'}) ...
                            || (strcmp(flag,'all') && D.zslice.header(chgdim).n==prevsz(chgdim))
                        flag = 'chgdata'; % no change in size
                    else
                        flag = 'chgdata&blocksize';
                    end
                    updateDisplay(D,flag)
                elseif ~chgclip
                    % the grid arrangement changes
                    switch flag
                        case 'chg'  % check this case first because it is this one that occurs when going fast through a list
                            updateDisplay(D,'chg',chgdim,e.ind)
                        case {'new' 'remove' 'perm'}
                            updateDisplay(D,flag,chgdim,e.ind)
                        case 'all'
                            ncur = size(D.htransform,chgdim);
                            n = D.zslice.sz(chgdim);
                            if n==ncur
                                updateDisplay(D,'chgdata')
                            elseif n>ncur
                                updateDisplay(D,'new',chgdim,ncur+1:n)
                                updateDisplay(D,'chg',chgdim,1:ncur)
                            else
                                updateDisplay(D,'remove',chgdim,n+1:ncur)
                                updateDisplay(D,'chg',chgdim,1:n)
                            end
                        case 'chg&new'
                            updateDisplay(D,'new',chgdim,e.ind{2})
                            updateDisplay(D,'chg',chgdim,e.ind{1})
                        case 'chg&rm'
                            updateDisplay(D,'remove',chgdim,e.ind{2})
                            updateDisplay(D,'chg',chgdim,e.ind{1})
                        otherwise
                            error('flag ''%s'' is not handled','flag')
                    end
                elseif chgclip
                    % all grid elements need to be updated
                    ncur = size(D.htransform,chgdim);
                    n = D.zslice.sz(chgdim);
                    if n==ncur
                        updateDisplay(D,'chgdata')
                    elseif n>ncur
                        updateDisplay(D,'chgdata')
                        updateDisplay(D,'new',chgdim,ncur+1:n)
                    else
                        updateDisplay(D,'remove',chgdim,n+1:ncur)
                        updateDisplay(D,'chgdata')
                    end
                end
            end
                      
           % reposition cross
           D.navigation.repositionCross()
           
           % update selections display
           D.navigation.displayselection('referentialchanged')

            % Update legend
            if strcmp(flag,'global')
                D.colordimID = [];
            elseif ~strcmp(flag,'chgdata') && isequal(chgdim,D.colordim)
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
                D.graph.setValueTicks()
                D.labels.updateLabels()                
                updateDisplay(D,'pos')
                D.navigation.repositionCross()
                D.navigation.displayselection()
            end
        end
        function axisresize(D)
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            D.graph.computeSteps()
            D.graph.setTicks()
            D.graph.setValueTicks()
            updateDisplay(D,'pos')
            D.navigation.repositionCross()
            D.navigation.displayselection()
        end
        function set.displaymode(D,mode)
            c = disableListener(D.listeners.axsiz); %#ok<MCSUP,NASGU> % prevent display update following automatic change of axis position
            % set property
            D.displaymode = mode;
            % for 'image' mode, check that layout is valid, and modify it if
            % necessary
            if strcmp(mode,'image') && length(D.layoutID.mergeddata)>1
                % it is not possible to have more than 1 color dimension ->
                % move other dimensions at 'mergeddata' location to 'y'
                newlayoutID = D.layoutIDmem;
                dimIDmove = D.layoutID.mergeddata(2:end);
                newlayoutID.mergeddata = setdiff(newlayoutID.mergeddata, dimIDmove, 'stable');
                newlayoutID.y = [newlayoutID.y dimIDmove];
                D.setLayoutID(newlayoutID) % this automatically updates display among other things
            else
                % update display
                D.checkzslicesize() % is zslice too large for being displayed
                D.graph.computeSteps() %#ok<MCSUP>
                D.graph.setTicks() %#ok<MCSUP>
                D.labels.updateLabels() %#ok<MCSUP>
                updateDisplay(D) % will call setValueTicks
                D.navigation.displayselection() %#ok<MCSUP>
            end
            % show/hide color legend
            displayColorLegend(D)
        end
    end
    
end