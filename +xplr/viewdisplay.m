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
        grid  % containers for line/images that will be translated/scaled
        hdisplay    % handles of line/images
        gridclip    % clipping for each grid cell
        signals_baseline    % subtracted baseline value for every signal if any
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
        external_dim
        external_dimID
        external_size
        overlap_dim
        clip_dim
        current_clip
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
            D.clipping = D.addComponent(xplr.cliptool(D)); % creates a menu
            
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
            if isempty(org.x)
                dim = [];
            else
                dim = org.x(1);
            end
            if strcmp(D.displaymode,'image')
                if ~isempty(org.y)
                    dim = [dim org.y(1) org.mergeddata];
                elseif ~isempty(org.mergeddata)
                    dim = [dim org.mergeddata];
                end
            end
            dim = sort(dim);
        end
        function dimID = get.internal_dimID(D)
            dim = D.internal_dim;
            dimID = [D.slice.header(dim).dimID];
        end
        function dim = get.external_dim(D)
            org = D.layout;
            dim = [org.x(2:end) org.y(1+strcmp(D.displaymode,'image'):end) org.xy org.yx];
            dim = sort(dim);
        end
        function dimID = get.external_dimID(D)
            dim = D.external_dim;
            dimID = [D.slice.header(dim).dimID];
        end
        function sz = get.external_size(D)
            sz = [D.slice.header(D.external_dim).n];
        end
        function dim = get.overlap_dim(D)
            % 'overlap' dimensions are those where several time course
            % signals can appear superimposed within a single grid cell
            if strcmp(D.displaymode,'time courses')
                dim = D.layout.mergeddata;
            else
                dim = zeros(1,0);
            end
        end
        function dim = get.clip_dim(D)
            % mergeddata and external dimensions together: these are the
            % dimensions where clipping range can be defined independently
            org = D.layout;
            dim = [org.x(2:end) org.y(1+strcmp(D.displaymode,'image'):end) org.xy org.yx org.mergeddata];
            dim = sort(dim);
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
            
            % when moving dimension label, immediate display update?
            fn_propcontrol(D.labels,'doImmediateDisplay','menu', ...
                {'parent',m,'label','Immediate display update when moving dimensions','separator','on'});
            
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
            if isequal(newlayoutID,D.layoutIDall), return, end
            c = disableListener(D.listeners.axsiz); %#ok<NASGU> % prevent display update following automatic change of axis position
            if nargin<3
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
        function clip = get.current_clip(D)
            % get clipping range of the current grid cell
            ijk = D.navigation.getPointIndexPosition('clip','round');
            ijk = row(round(D.graph.slice2zslice(ijk)));
            clip = subsref_dim(D.gridclip,1+D.clip_dim,ijk(D.clip_dim));
            % center if 'adjust by the mean' mode
            if strcmp(D.displaymode,'time courses') && ~isempty(D.clipping.align_signals)
                clip = clip - feval(D.clipping.align_signals,clip,1);
            end
        end
        function setClip(D,clip,all_cells)
            % input
            if ~isnumeric(clip) || length(clip)~=2 || diff(clip)<=0 || any(isnan(clip)|isinf(clip))
                xplr.debuginfo('stop','clip value is not valid')
                return
            end
            if nargin<3, all_cells = false; end
            % set property
            if all_cells
                D.gridclip(1,:) = clip(1);
                D.gridclip(2,:) = clip(2);
            else
                ijk = D.navigation.getPointIndexPosition('clip','round');
                ijk = row(round(D.graph.slice2zslice(ijk)));
                dim_indp = D.clipping.independent_dim;
                D.gridclip = subsasgn_dim(D.gridclip,[1 1+dim_indp],[1 ijk(dim_indp)],clip(1));
                D.gridclip = subsasgn_dim(D.gridclip,[1 1+dim_indp],[2 ijk(dim_indp)],clip(2));
            end
            % update display; note that this might also modify D.gridclip
            % if D.cliping.align_signals is not empty
            updateDisplay(D,'clip')
        end
        function autoClip(D, all_cells)
            if nargin<2, all_cells = false; end
            if all_cells
                D.gridclip(:) = NaN;
            else
                ijk = D.navigation.getPointIndexPosition('clip','round');
                ijk = row(round(D.graph.slice2zslice(ijk)));
                dim_indp = D.clipping.independent_dim;
                D.gridclip = subsasgn_dim(D.gridclip,1+dim_indp,ijk(dim_indp),NaN);
            end
            updateDisplay(D,'clip')            
        end
        function clip = get_clip_range(D,data)
            clip = fn_clip(data(:),D.clipping.autoclipmode,'getrange');
            if isinf(clip(1)), clip(1) = -1e6; end
            if isinf(clip(2)), clip(2) = 1e6; end
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
            %---
            % available flags are: 'clip' 'global' 'chgdata'
            % 'chgdata&blocksize' 'new' 'remove' 'chg' 'perm' 'pos' 'color' 
            if nargin<3, dim = []; end
            
            % Is data too large for being displayed?
            if D.nodisplay
                % too many grid elements: cancel display!
                deleteValid(D.grid) % this will also delete children D.hdisplay
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
            sz = D.zslice.sz;
            org = D.layout; % convert from dim ID to dim numbers
            internaldim = D.internal_dim;
            externaldim = D.external_dim;
            overlapdim = D.overlap_dim;
            clipdim = D.clip_dim;
            elementsdim = sort([overlapdim externaldim]);
            
            % What to do
            if nargin<2, flag = 'global'; end
            if ~fn_ismemberstr(flag,{'clip' 'global' 'chgdata' 'chgdata&blocksize' 'new' 'remove' 'chg' 'perm' 'pos' 'color'})
                error 'flag not handled'
            end
            doreset = strcmp(flag,'global');
            dim_external = ismember(dim,externaldim);
            donew = strcmp(flag,'new');
            doremove = strcmp(flag,'remove');
            doposition = ~fn_ismemberstr(flag,{'chg' 'chgdata' 'clip' 'color'});
            dodataall = fn_ismemberstr(flag,{'clip' 'global' 'chgdata' 'chgdata&blocksize' 'perm' 'color'}); % color is set when updating ydata, but updating ydata is actually not necessary when only color changes...
            dodataselect = fn_ismemberstr(flag,{'new' 'chg'});
            dodata = dodataall || dodataselect;
            dodataselect_grid = dodataselect && dim_external;
            dodataselect_overlap = dodataselect && ~dim_external;
            if dodataselect, assert(dodataselect_grid || dodataselect_overlap), end
            dochgx = strcmp(flag,'chgdata&blocksize');
            docolor = dotimecourses && ~fn_ismemberstr(flag,{'chg' 'chgdata' 'chgdata&blocksize' 'clip'});
            
            % Grid size
            grid_size = ones(1,max(D.nd,2)); 
            grid_size(externaldim) = sz(externaldim);
            hdisplay_size = ones(1,max(D.nd,2)); 
            hdisplay_size(elementsdim) = sz(elementsdim);
            gridclip_size = ones(1,max(D.nd,2)); 
            gridclip_size(clipdim) = sz(clipdim);
            
            % Check that current grid are valid
            prev_size_check = grid_size;
            if ~isempty(dim), prev_size_check(dim) = prev_size_check(dim)+(doremove-donew)*length(ind); end
            if ~isequal(strictsize(D.grid,length(grid_size)),prev_size_check) ...
                    || ~all(ishandle(D.grid(:)))
                [doreset, doposition, dodataall] = deal(true);
                docolor = dotimecourses;
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
            if doreset          % reset display and grid elements
                deleteValid(D.grid) % this will also delete children D.hdisplay
                D.grid = gobjects(grid_size);
                D.gridclip = NaN([2 gridclip_size]);
                D.hdisplay = gobjects(hdisplay_size);
                [doposition, dodataall] = deal(true);
            elseif donew      	% new grid elements
                if ismember(dim,externaldim)
                    D.grid = subsasgn_dim(D.grid,dim,ind,gobjects);
                end
                D.gridclip = subsasgn_dim(D.gridclip,1+dim,ind,NaN);
                D.hdisplay = subsasgn_dim(D.hdisplay,dim,ind,gobjects);
            elseif doremove     % remove grid elements
                if ismember(dim,externaldim)
                    deleteValid(subsref_dim(D.grid,dim,ind)) % this also deletes the children hdisplay objects
                    D.grid = subsasgn_dim(D.grid,dim,ind,[]);
                end
                D.gridclip = subsasgn_dim(D.gridclip,1+dim,ind,[]);
                deleteValid(subsref_dim(D.hdisplay,dim,ind)) % this also deletes the children hdisplay objects
                D.hdisplay = subsasgn_dim(D.hdisplay,dim,ind,[]);
            elseif strcmp(flag,'chgdata&blocksize') && strcmp(D.displaymode,'image') && ismember(dim,org.mergeddata)
                % D.gridclip might need to be changed in dimension
                % org.mergeddata for images
                % TODO: we have too many cases with this mergeddata, need
                % to simplify somehow...
                n_prev = min(size(D.gridclip,1+dim),4);
                n_cur = min(sz(dim),4);
                if n_cur > n_prev
                    D.gridclip = subsasgn_dim(D.gridclip, 1+dim, n_prev+1:n_cur, NaN);
                elseif n_cur < n_prev
                    D.gridclip = subsasgn_dim(D.gridclip, 1+dim, n_cur+1:n_prev, []);
                end
            end
            
            % Update clipping: if clipping is independent for each element,
            % postpone clip computation to when these elements will be
            % extracted and displayed. Otherwise compute now.
            if dodata
                displayed_data = D.zslice.data;
                % correct data to align signals on their mean or median?
                if strcmp(D.displaymode,'time courses')
                    align_signals = D.clipping.align_signals;
                else
                    align_signals = false;
                end
                if align_signals
                    D.signals_baseline = displayed_data;
                    if ~isempty(org.x)
                        D.signals_baseline = feval(align_signals,D.signals_baseline,org.x(1)); 
                        displayed_data = fn_subtract(displayed_data, D.signals_baseline);
                    end
                else
                    D.signals_baseline = [];
                end
                % invalidate clip values that need to be recomputed
                if D.clipping.adjust_to_view && ~strcmp(flag,'clip')
                    if ismember(dim,D.clipping.independent_dim)
                        D.gridclip = subsasgn_dim(D.gridclip,1+dim,ind,NaN);
                    else
                        D.gridclip(:) = NaN;
                    end
                end
                % compute clip
                dim_indpc = D.clipping.independent_dim; % dimensions with independent clipping
                clip_at_elements_level = isequal(elementsdim, dim_indpc);
                if strcmp(flag,'perm')
                    % mere permutation
                    D.gridclip = subsref_dim(D.gridclip,dim,ind);
                elseif clip_at_elements_level
                    % every element will be clipped independently: postpone
                    % clip computation to when elements will be extracted,
                    % to avoid unneeded repetitive extractions
                else
                    % compute clip independently in sub-parts of the data
                    sz_indpc = ones(1,D.nd); sz_indpc(dim_indpc) = sz(dim_indpc);
                    for idx_indpc = 1:prod(sz_indpc)
                        ijk_indpc = row(fn_indices(sz_indpc(dim_indpc),idx_indpc,'g2i'));
                        dat_part = subsref_dim(displayed_data,dim_indpc,ijk_indpc);
                        clip_part = subsref_dim(D.gridclip,1+dim_indpc,ijk_indpc);
                        idx = find(~isnan(clip_part(1,:)),1);
                        if isempty(idx)
                            % need to recompute
                            clip = D.get_clip_range(dat_part);
                        else
                            % use current clip value where it is not
                            % defined
                            clip = clip_part(:,idx);
                        end
                        D.gridclip = subsasgn_dim(D.gridclip,[1 1+dim_indpc],[1 ijk_indpc],clip(1));
                        D.gridclip = subsasgn_dim(D.gridclip,[1 1+dim_indpc],[2 ijk_indpc],clip(2));
                    end
                end
            end
            
            % Create or modify grid containers
            if doposition
                % (list of grid cell indices)
                idx_grid_list = 1:prod(grid_size);
                if dodataselect_grid && ~doposition
                    % not all grid elements need to be visited ('chg' flag)
                    idx_grid_list = reshape(idx_grid_list,grid_size);
                    subs = substruct('()',repmat({':'},1,D.zslice.nd));
                    subs.subs{dim} = ind;
                    idx_grid_list = row(subsref(idx_grid_list,subs));
                end
                ijk_grid_list = fn_indices(grid_size,idx_grid_list,'g2i');
                % (prepare dispatch)
                M = D.graph.gettransform(ijk_grid_list);
                % (loop on grid cells)
                for u = 1:length(idx_grid_list)
                    idx_grid = idx_grid_list(u);
                    ijk_grid = ijk_grid_list(:,u);
                    dodata_cell = dodata && (~dodataselect_grid || any(ind==ijk_grid(dim)));
                    docreatecur = doreset || (donew && dodata_cell);
                    if docreatecur
                        D.grid(idx_grid) = hgtransform('parent',D.ha,'matrix',M(:,:,idx_grid),'HitTest','off');
                    else
                        set(D.grid(idx_grid),'matrix',M(:,:,idx_grid))
                    end
                end
            end
            
            % Create or modify individual signals/images
            if dodata
                % list of all signals/images indices
                idx_hdisplay_list = reshape(1:prod(hdisplay_size),hdisplay_size);
                if dodataselect
                    % not all objects need to be visited
                    subs = substruct('()',repmat({':'},1,D.nd));
                    subs.subs{dim} = ind;
                    idx_hdisplay_list = subsref(idx_hdisplay_list,subs);
                end
                % (reorganize list as: rows=different grid cells, columns=overlapped elements inside same grid cell)
                idx_hdisplay_list = fn_reshapepermute(idx_hdisplay_list,{externaldim overlapdim internaldim});
                [ngrid, noverlap] = size(idx_hdisplay_list); % can be smaller than total numbers of grids/overlaps
                ijk_hdisplay_list = fn_indices(hdisplay_size,idx_hdisplay_list(:),'g2i');
                ijk_hdisplay_list = reshape(ijk_hdisplay_list,[D.nd ngrid noverlap]);
                % (corresponding grid cells) 
                ijk_grid_list = ijk_hdisplay_list(:,:,1);
                ijk_grid_list(overlapdim,:) = 1;
                idx_grid_list = fn_indices(grid_size,ijk_grid_list,'i2g');
                
                % subs structure for slicing
                subs = substruct('()',repmat({':'},1,length(sz)));
                % (and permutation that follows)
                internalperm = zeros(1,D.nd);
                if strcmp(D.displaymode,'image')
                    % reorder dimensions as y-x-channel-others
                    if ~isempty(org.y), internalperm(1) = org.y(1); end
                    if ~isempty(org.x), internalperm(2) = org.x(1); end
                    if ~isempty(org.mergeddata), internalperm(3) = org.mergeddata; end
                else
                    % reorder dimensions as x-others
                    if ~isempty(org.x), internalperm(1) = org.x(1); end
                end
                internalperm(~internalperm) = setdiff(1:D.nd,internaldim);
                
                % go! loop on grid cells and on elements inside these cells
                for u = 1:ngrid
                    idx_grid = idx_grid_list(u);
                    for v = 1:noverlap
                        idx_hdisplay = idx_hdisplay_list(u,v);
                        ijk_hdisplay = ijk_hdisplay_list(:,u,v);

                        % get the data and compute clipping if necessary
                        [subs.subs{elementsdim}] = dealc(ijk_hdisplay(elementsdim));
                        xi = subsref(displayed_data, subs);
                        xi = permute(xi, internalperm);
                        xi = fn_float(xi);
                        clipi = subsref_dim(D.gridclip,1+elementsdim,ijk_hdisplay(elementsdim));
                        if clip_at_elements_level && any(isnan(clipi))
                            % independent clip for this element
                            clipi = D.get_clip_range(xi);
                            clipi = repmat(clipi(:),[1 size(xi,3)]); 
                            D.gridclip = subsasgn_dim(D.gridclip,1+elementsdim,ijk_hdisplay(elementsdim),clipi);
                        end
                        % display it
                        if dotimecourses
                            xi = (xi-clipi(1))/diff(clipi);
                            nt = size(xi,1);
                            if nt==1
                                lineopt = {'linestyle','none','marker','.'};
                            else
                                lineopt = {'linestyle','-','marker','none'};
                            end
                            if ~ishandle(D.hdisplay(idx_hdisplay))
                                hl = line(1:nt,xi, ...
                                    'parent',D.grid(idx_grid),'HitTest','off',lineopt{:});
                                D.hdisplay(idx_hdisplay) = hl;
                            else
                                hl = D.hdisplay(idx_hdisplay);
                                if dochgx
                                    set(hl,'xdata',1:nt,lineopt{:})
                                end
                                set(hl,'ydata',xi)
                            end
                            if docolor, set(hl,'color',cmap(ijk_hdisplay(cdim),:)), end
                        else
                            % size in color dimension must be 1, 3 or 4;
                            % correct if it is not the case
                            alpha = [];
                            nc = size(xi,3);
                            if nc > 4
                                xi = xi(:,:,1:4);
                            end
                            im = D.colormap.color_image(xi, clipi);
                            if nc == 2
                                % add a third blue channel set to zero
                                im(:,:,3) = 0;
                            elseif nc == 4
                                [im, alpha] = deal(im(:,:,1:3), im(:,:,4));
                            end
                            if ~ishandle(D.hdisplay(idx_hdisplay))
                                % y coordinates are negative to orient the
                                % image downward (see also comment inside of
                                % displaygaph.gettransform method, where the
                                % y-scale of the hgtransform cannot be set
                                % negative)
                                D.hdisplay(idx_hdisplay) = surface([.5 size(im,2)+.5],[-.5 -.5-size(im,1)],zeros(2), ...
                                    'parent',D.grid(idx_grid), ...
                                    'EdgeColor','none','FaceColor','texturemap', ...
                                    'CDataMapping','scaled','FaceAlpha','texturemap', ...
                                    'CData',im,'AlphaData',alpha, ...
                                    'HitTest','off');
                            elseif dochgx
                                set(D.hdisplay(idx_hdisplay),'xdata',[.5 size(im,2)+.5],'ydata',[-.5 -.5-size(im,1)], ...
                                    'CData',im,'AlphaData',alpha)
                            else
                                set(D.hdisplay(idx_hdisplay),'CData',im,'AlphaData',alpha)
                            end
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
                uistack(D.grid(:), 'bottom')
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
            
            % Update display
            if fn_ismemberstr(flag,{'global' 'chgdim'})
                % Reset display
                updateDisplay(D,'global')
            elseif strcmp(flag,'chgdata')
                % No change in size, all data need to be redisplayed
                updateDisplay(D,'chgdata')
            else
                % Smart display update
                if ismember(chgdimID, D.internal_dimID)
                    % changes are within elements (the grid arrangement
                    % remains the same)
                    if fn_ismemberstr(flag,{'perm' 'chg'}) ...
                            || (strcmp(flag,'all') && D.zslice.header(chgdim).n==prevsz(chgdim))
                        flag = 'chgdata'; % no change in size
                    else
                        flag = 'chgdata&blocksize';
                    end
                    updateDisplay(D,flag,chgdim,e.ind)
                else
                    % the grid arrangement changes
                    switch flag
                        case {'chg' 'new' 'remove' 'perm'}
                            updateDisplay(D,flag,chgdim,e.ind)
                        case 'all'
                            ncur = size(D.grid,chgdim);
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
                            updateDisplay(D,'chg',chgdim,e.ind{1})
                            updateDisplay(D,'new',chgdim,e.ind{2})
                        case 'chg&rm'
                            updateDisplay(D,'chg',chgdim,e.ind{1})
                            updateDisplay(D,'remove',chgdim,e.ind{2})
                        otherwise
                            error('flag ''%s'' is not handled','flag')
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