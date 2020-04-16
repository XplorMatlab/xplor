classdef displaynavigation < xplr.graphnode
% display navigation

    properties (SetAccess='private')
        D                                       % parent xplr.viewdisplay
        ha = gobjects
        hf = gobjects
        graph
        crossCenter
        crossDataValue                          % data value at cross position
        cross = gobjects                        % display cross selector
        sliders = struct('x',[],'y',[]);        % slider objects
        zoomfilters = struct('x',[],'y',[]);    % connected zoom filters
        pointfilters = {};
    end
    properties
        selectionfilter         % filter being modified by the selections
        selectiondisplay        % displays of selectionnd: for each selection, 2 lines and a text
        selectionsavefile       % name of file for saving current selection
        selection_context       
    end
    properties (SetObservable, AbortSet=true)
        showcross = true;
        crosscolor = [0 0 0];
        crossalpha = .5;
        selectiondimID          % dimensions to which these selections apply, identified by its id
        selection2Dshape = 'ellipse'; % 'poly', 'free', 'rect', 'ellipse', 'ring', 'segment', 'openpoly', 'freeline'
        selectionround1Dmeasure = true; 
        selectionshow = 'shape+name';
        selectionedit = false
        selectionpromptname = false
    end
    properties (Dependent)
        selectiondim            % dimension(s) to which selections apply, identified by its(their) number
    end
    properties (Dependent, SetAccess='private')
        selection               % list of selectionnd object from the filter
    end

    % Constructor
    methods
        function N = displaynavigation(D)
            % parent xplr.viewdisplay object and other external objects
            N.D = D;
            N.ha = D.ha;
            N.hf = D.V.hf;
            N.graph = D.graph;
            
            % buttons
            init_buttons(N)
            
            % cross
            N.displaycross()
            
            % sliders
            init_sliders(N)
            
            % data value display
            init_value_display(N)

            % connect sliders to the active dimensions of the display
            % (note that this is in fact redundant with call in
            % viewdisplay.slicechange when viewdisplay object is created)
            connectZoomFilter(N)

            % axes_click actions
            set(D.ha,'buttondownfcn',@(u,e)axes_click(N))

            % scroll wheel zooming
            fn_scrollwheelregister(D.ha,@(n)N.Scroll(n))
            
            % selection menu
            uimenu(N.hf,'Label','Selection','callback',@(m,e)N.selectionMenu(m))
        end
        function init_buttons(N)
        % function init_buttons(N)
        % 3 buttons that control clipping
            
            % first button to adjust clipping with axes_click movements:
            % display image on it indicating how image luminance and
            % contrast change upon axes_click movements
            [ii jj] = meshgrid(-13:0,13:-1:0); x=(0-ii)./(jj-ii)-.5; x(end)=0;
            u = uicontrol('parent',N.D.hp, ...
                'enable','inactive','cdata',fn_clip(sin(pi*x),[-1 1],'gray'), ...
                'buttondownfcn',@(u,e)moveclip(N));
            fn_controlpositions(u,N.ha,[1 1 0 0],[-1 -16 16 16])

            % two next buttons control extent of clipping
            u = uicontrol('parent',N.D.hp, ...
                'string','+','fontsize',8, ...
                'callback',@(u,e)cliprange(N,'+'));
            fn_controlpositions(u,N.ha,[1 1 0 0],[-1 -32 16 16])
            u = uicontrol('parent',N.D.hp, ...
                'string','-','fontsize',8, ...
                'callback',@(u,e)cliprange(N,'-'));
            fn_controlpositions(u,N.ha,[1 1 0 0],[-1 -48 16 16])
        end
        function init_sliders(N)
            N.sliders.x = fn_slider('parent',N.D.hp,'mode','area', ...
                'layout','right','callback',@(u,evnt)chgzoom(N,'x',u));
            N.sliders.y = fn_slider('parent',N.D.hp,'mode','area', ...
                'layout','down','callback',@(u,evnt)chgzoom(N,'y',u));
            pcol = get(N.D.hp,'backgroundcolor');
            set([N.sliders.x N.sliders.y],'visible','off','scrollwheel','on','value',[0 1], ...
                'backgroundcolor',pcol*.75,'slidercolor',pcol*.95)
            fn_controlpositions(N.sliders.x,N.ha,[0 1 1 0], [0 0 0 12]);
            fn_controlpositions(N.sliders.y,N.ha,[1 0 0 1], [0 0 12 -48]);
        end
        function init_value_display(N)
            N.crossDataValue = uicontrol('Parent',N.D.hp,'style','text','enable','inactive', ...
                'fontsize',8,'horizontalalignment','right');
            fn_controlpositions(N.crossDataValue,N.D.hp,[1 0],[-75 10 65 15])
        end
    end
    
    % Get/Set dependent
    methods
        function d = get.selectiondim(N)
            d = N.D.slice.dimensionNumber(N.selectiondimID);
        end
        function set.selectiondim(N,dim)
            N.selectiondimID = N.D.slice.dimensionID(dim);
        end
    end
    
    % Clip control
    methods
        function moveclip(N)
            switch get(N.hf,'selectiontype')
                case 'normal'       % change clip
                    clip0 = N.D.clip;
                    e0 = diff(clip0);
                    clipcenter = N.D.clipping.center;
                    switch N.D.clipping.adjust
                        case 'none'
                            % nothing
                        otherwise
                            clipcenter = 0;
                    end
                    if ~isempty(clipcenter), clip0 = clipcenter + [-.5 .5]*e0; end
                    p0 = get(N.hf,'currentpoint');
                    ht = uicontrol('style','text','position',[2 2 200 17],'parent',N.hf);
                    % change clip
                    moveclipsub() % this displays the bottom-left numbers
                    fn_buttonmotion(@moveclipsub,N.hf)
                    delete(ht)
                case 'open'         % use default clipping
                    autoClip(N.D)
            end
            function moveclipsub
                % 'naive' new clip
                p = get(N.hf,'currentpoint');
                dp = p-p0;
                if ~isempty(clipcenter), dp = [-1 1]*(dp(2)-dp(1))/2; end
                FACT = 1/100;
                clip = clip0 + dp*(e0*FACT);
                % it might be that we have diff(clip)<=0 here! apply some
                % transformation to solve that
                e = diff(clip);
                thr = e0/10;
                if e<thr
                    %e = thr*exp(e/thr-1); % goes from thr for e=thr to 0 for e=-Inf
                    e = thr^2 / (2*thr-e); % goes from thr for e=thr to 0 for e=-Inf
                    clip = mean(clip) + [-.5 .5]*e;
                end     
                % update display
                set(ht,'string',sprintf('min: %.3f,  max: %.3f',clip(1),clip(2)))
                N.D.setClip(clip)
            end
        end
        function cliprange(N,flag)
            % current clip extent
            clip = N.D.clip;
            m = mean(clip);
            e = diff(clip);
            
            % round it to a nice value
            e10 = 10^floor(log10(e));
            e = e / e10;
            vals = [.75 1 1.5 2 3 4 5 7.5 10 15];
            f = find(e*1.1>vals,1,'last');
            
            % update as specified
            f = f + fn_switch(flag,'+',-1,'-',1);
            e = e10 * vals(f);
            clip = m + [-.5 .5]*e;
            
            % set clip
            N.D.setClip(clip)
        end
    end
    
    % axes_click actions
    methods
        function axes_click(N, flag)
            % function axes_click(N)
            % function axes_click(N, 'pointonly')
            %---
            % Options:
            % - 'pointonly' set this flag when we have already clicked
            %               and released the axes_click button (for example
            %               because we clicked on the cross): in this case
            %               do not start drag zooming or ROI selection
            
            if nargin<2, flag = ''; end
            pointonly = strcmp(flag,'pointonly');
            point =  get(N.D.ha,'CurrentPoint'); point = point(1,[1 2])';
            activedim = [N.D.activedim.x N.D.activedim.y];
            switch get(N.hf,'SelectionType')
                case 'normal'
                    % zoom in or select point
                    if pointonly                        
                        dozoom = false;
                    else
                        rect = fn_mouse(N.ha,'rectaxp-');
                        dozoom = any(any(abs(diff(rect,1,2))>1e-2));
                    end
                    if dozoom
                        if isempty(activedim)
                            % determine in which dimension to zoom as the most
                            % exterior dimension(s) where at least 2
                            % elements are selected
                            zijk = N.graph.graph2zslice(rect);
                            ijk = N.graph.graph2slice(rect);
                            nonsingleton = logical(diff(round(zijk),1,2)); % nd*1 vector
                            org = N.D.layout;
                            xydim = [org.xy org.yx];
                            if nonsingleton(xydim)
                                zoomdim = xydim;
                            else
                                xdim = org.x(find(nonsingleton(org.x),1,'last'));
                                ydim = org.y(find(nonsingleton(org.y),1,'last'));
                                zoomdim = [xdim ydim];
                            end
                            zoom = ijk(zoomdim,:)';
                        else
                            zoomdim = activedim;
                            zoom = N.graph.graph2slice(rect,'subdim',zoomdim)';
                        end
                        for i=1:length(zoomdim), zoom(:,i) = sort(zoom(:,i)); end
                        N.D.zoomslicer.setZoom(zoomdim,zoom)
                    else
                        N.manualclickmovecross(point);
                    end
                case 'open'
                    % zoom reset
                    if isempty(activedim)
                        zoomoutdim = 1:N.D.nd;
                    else
                        zoomoutdim = activedim;
                    end
                    zoom = repmat(':',1,length(zoomoutdim));
                    N.D.zoomslicer.setZoom(zoomoutdim,zoom)
                case 'alt'
                    % on right click: create a new selection in N.selection
                    % depending on the parameter selectionshape
                    if ~isempty(N.selectiondimID)
                        selslice = N.selectionMouse();
                        if isempty(selslice), return, end

                        % prompt for selection name
                        options = {};
                        if N.selectionpromptname
                            name = inputdlg('Selection name','xplor');
                            if ~isempty(name), options = {'Name' name{1}}; end
                        end
                        
                        % update filter
                        N.selectionfilter.updateSelection('new',selslice,options{:})
                    end
            end
        end
    end
    
    % Cross point selection
    methods
        function connectPointFilter(N,dim,key)
            if nargin < 3, key = 1; end

            % Disconnect current point filters
            if nargin < 2
                disconnectPointFilter(N)
                dim = 1:N.D.slice.nd;
            else
                disconnectPointFilter(N,dim)
            end
            
            % Connect new point filters: loop on dimensions
            for d = dim
                linkkey = key;
                head = N.D.slice.header(d);
                % no interest in creating and controling a filter for a
                % dimension with only 1 value
                if head.n ==1 
                    N.pointfilters{d} = [];
                    continue
                end
                % get filter from bank or create one for the header in this
                % dimension
                P = xplr.bank.getPointFilter(linkkey,head,N); % FilterAndPoint filter
                N.pointfilters{d} = P;
                % listen to the point filter event
                N.addListener(P,'ChangedOperation',@(u,e)movedPoint(N,e))
            end
        end
        function disconnectPointFilter(N,dim)
            if nargin < 2
                dim = 1:length(N.pointfilters);
            end
            for d = dim
                P = N.pointfilters{d};
                if isempty(P), continue, end
                xplr.bank.unregisterFilter(P,N)
                N.disconnect(P)  % this is the same as F.disconnect(N)!
            end
        end
        function movedPoint(N,e)
            if ~strcmp(e.type,'point'), return, end
            ijk = getPointIndexPosition(N);
            N.crossCenter = N.graph.slice2graph(ijk);
            update_cross_visibility(N);
        end
        function ijk = getPointIndexPosition(N)
            nd = N.D.slice.nd;
            ijk = ones(nd, 1);
            for d = 1:nd
                P = N.pointfilters{d};
                if isempty(P), continue, end
                ijk(d) = P.index0;
            end
        end
        function repositionCross(N)
            ijk = getPointIndexPosition(N);
            N.crossCenter = N.graph.slice2graph(ijk);
            update_cross_visibility(N);
        end
        function displaycross(N)
           
            % cross
            N.cross(1) = line('Parent',N.D.ha,'ydata',[-.5 .5]);
            N.cross(2) = line('Parent',N.D.ha,'xdata',[-.5 .5]);
            N.cross(3) = line('Parent',N.D.ha,'xdata',[0 0],'ydata',[0 0]); % a single point
            set(N.cross,'Color',[N.crosscolor N.crossalpha]) % cross is semi-transparent!
            
            % position
            N.crossCenter = [0 0];
            
            % callbacks
            for i=1:3
                set(N.cross(i),'buttondownfcn',@(u,e)manualmovecross(N,i))
            end
        end
        function set.crossCenter(N, crossCenter)

            % re-display cross if it was deleted (happens upon
            % D.resetDisplay)
            if ~all(isvalid(N.cross))
                deleteValid(N.cross)
                N.displaycross()
            end
            
            % set property
            N.crossCenter = crossCenter;

            % move the cross
            set(N.cross(1),'XData',crossCenter([1 1]))
            set(N.cross(2),'YData',crossCenter([2 2]))
            set(N.cross(3),'XData',crossCenter([1 1]),'YData',crossCenter([2 2]))

        end
        function manualmovecross(N,il)
            if ~ismember(get(N.hf,'selectiontype'),{'normal' 'open'})
                % not a left click: execute callback for axes
                axes_click(N)
                return
            end
            set(N.hf,'pointer',fn_switch(il,1,'left',2,'top',3,'cross'))
            
            % prepare a time to pan zoom while moving the cross!
            do_drag_zoom = (il == 1) && isequal(N.D.layout.x,1) && isequal(N.D.activedim.x,1);
            if do_drag_zoom
                Z = N.D.zoomfilters(1);
                do_drag_zoom = ~strcmp(Z.zoom,':');
            end
            if do_drag_zoom
                drag_timer = timer('timerfcn',@drag_zoom,'ExecutionMode','fixedSpacing','period',.01);
                zoom_width = diff(Z.zoom);
                zoom_min = .5;
                zoom_max = .5+Z.headerin.n;
                drag_zone = .15; % width of special zone where zoom dragging occurs
                drag_speed = 0;
            end
            point = [];
            
            anymove = fn_buttonmotion(@movecrosssub,N.hf,'moved?');            
            set(N.hf,'pointer','arrow')
            if do_drag_zoom
                stop(drag_timer)
            end
            if ~anymove
                % execute callback for axes
                axes_click(N, 'pointonly')
                return
            end
            
            function movecrosssub
                %anymove = true;
                p = get(N.D.ha,'currentpoint'); p = p(1,1:2);
                switch il
                    case 1
                        point = [p(1) N.crossCenter(2)];
                    case 2
                        point = [N.crossCenter(1) p(2)];
                    case 3
                        point = p;
                end
                % move the cross
                manualclickmovecross(N,point)
                % pan view if we are close to edge!
                if do_drag_zoom
                    if abs(p(1))>.5-drag_zone
                        drag_speed = sign(p(1)) * ((abs(p(1)) - (.5-drag_zone))/drag_zone)/10;
                        if ~boolean(drag_timer.Running)
                            start(drag_timer)
                        end
                    else
                        stop(drag_timer)
                    end
                end
            end
            
            function drag_zoom(~,~)
                zoom = Z.zoom + drag_speed * zoom_width;
                if zoom(1) < zoom_min
                    zoom = zoom_min + [0 zoom_width];
                elseif zoom(2) > zoom_max
                    zoom = zoom_max - [zoom_width 0];
                end
                Z.setZoom(zoom)
                manualclickmovecross(N,point)
            end
            
        end        
        function manualclickmovecross(N,point)
            % move the cross to the selected point
            ijk = N.graph.graph2slice(point,'invertible',true);
            
            % round indices values in dimensions with categorical headers
            categorical = [N.D.slice.header.categorical];
            ijk(categorical) = round(ijk(categorical));
            
            % update the point filters (only for dimensions where the point
            % remains within the slice)
            for d = find(~isPointOutOfDisplay(N,point,true))
                P = N.pointfilters{d};
                if ~isempty(P)
                    P.index = ijk(d);
                end
            end
        end        
        function removecross(N)
            delete(N.cross)
        end      
        function update_cross_visibility(N)
            % do not show cross?
            if ~N.showcross
                set(N.cross,'Visible','off')
                return
            end           
                        
            % dims that are out of display
            ijk = getPointIndexPosition(N);
            dimOutOfDisplay = N.isIndexOutOfDisplay(ijk,true);
           
            % Hide the vertical bar if all dimensions on x are singletons or if
            % crossCenter is out of display on any dimension on x
            layout = N.D.layout;
            xdim = [layout.x layout.xy layout.yx];
            x_singleton = isempty(xdim);
            x_isOutOfDisplay = any(dimOutOfDisplay(xdim));
            N.cross(1).Visible = onoff(~x_singleton && ~x_isOutOfDisplay);
            
            % Same things for horizontal bar
            ydim = [layout.y layout.xy layout.yx];
            y_singleton = isempty(ydim);
            y_isOutOfDisplay = any(dimOutOfDisplay(ydim));
            N.cross(2).Visible = onoff(~(y_singleton|y_isOutOfDisplay));
            
            % Cross center
            updateCrossCenterVisibility(N);

            % Cross Value
            updateValueDisplay(N);
        end
        function updateCrossCenterVisibility(N)
            %  if one of the dimension of the cross is hidden, hide the
            % cross center as well
            N.cross(3).Visible = onoff(boolean(N.cross(1).Visible) && boolean(N.cross(2).Visible));
        end
        function updateValueDisplay(N)
            ijk = getPointIndexPosition(N);
            idx = fn_indices(N.D.slice.sz, round(ijk));
            value = N.D.slice.data(idx);

            % Test to display the value as "val(d1,d2,d3,...)=value"
            %set(N.crossDataValue,'String',['val(' num2str(ijk(1),'%.3g') ',' num2str(ijk(2),'%.3g') ')=' ...
            %            num2str(value,'%.3g')])

            set(N.crossDataValue,'String',num2str(value,'%.3g'))
            %disp("Cross value updated");
        end
        % cross color, transparency, and global visibility
        function set.showcross(N,value)
            N.showcross = value;
            N.update_cross_visibility()
        end
        function set.crosscolor(N,color)
            color = fn_colorbyname(color);
            if isempty(color), error 'wrong color', end
            N.crosscolor = color;
            set(N.cross,'color',[N.crosscolor N.crossalpha])
        end
        function set.crossalpha(N,alpha)
            N.crossalpha = alpha;
            set(N.cross,'color',[N.crosscolor N.crossalpha])
        end
    end
    
    % ROI selection
    methods
        function selectionMenu(N,m)
            % Selection menu is populated when being activated
            delete(get(m,'children'))
            
            header = N.D.slice.header;
            seldim = N.selectiondim;
            sellabels = {header(seldim).label};
            layout = N.D.layout;

            % Set selection dimension
            % (no active control? -> propose selection in the 'internal' dimensions)
            if isempty(N.selectiondimID) && ~isempty([layout.x layout.y])
                if ~isempty(layout.x) && (strcmp(N.D.displaymode,'time courses') || isempty(layout.y)) 
                    % time courses (or image without y dimension) -> 1D selection
                    propdim = layout.x(1);
                elseif strcmp(N.D.displaymode,'image') && isempty(layout.x)
                    % image but no x dimension -> 1D vertical selection
                    propdim = layout.y(1);
                else
                    % image -> 2D selection
                    propdim = [layout.x(1) layout.y(1)];
                end
                proplabels = {header(propdim).label};
                if isscalar(propdim)
                    itemlabel = ['Control selection in dimension ' proplabels{1}];
                else
                    itemlabel = ['Control selection in dimensions ' fn_strcat(proplabels,',')];
                end
                propdimID = [header(propdim).dimID];
                uimenu(m,'label',itemlabel, ...
                    'callback',@(u,e)set(N,'selectiondimID',propdimID))
            end
            % (sub-menu with other possible dimensions)
            switch length(N.selectiondimID)
                case 0
                    info = 'Control selection in dimension...';
                case 1
                    info = ['Control selection in dimension: ' sellabels{1}];
                case 2
                    info = ['Control selection in dimensions: ' fn_strcat(sellabels,',')];
                otherwise
                    error 'programming: selection control in more than 2 dimensions'
            end
            % (submenu for other possibilities)
            m2 = uimenu(m,'label',info);
            % (1D: dimension location must be x, y or xy)
            dimok = sort([layout.x layout.y layout.xy]);
            fn_propcontrol(N,'selectiondimID', ...
                {'menugroup' {header(dimok).dimID} {header(dimok).label}}, ...
                {'parent',m2});
            % (2D: dimension locations must be respectively x and y)
            available = cell(2, length(layout.x), length(layout.y));
            if ~isempty(available)
                % selections with first dim on x-axis, second dim on y-axis
                for i = 1:length(layout.x)
                    for j = 1:length(layout.y)
                        d = [layout.x(i) layout.y(j)];
                        available{1,i,j} = [header(d).dimID];
                        available{2,i,j} = fn_strcat({header(d).label},',');
                    end
                end
                fn_propcontrol(N,'selectiondimID', ...
                    {'menugroup' available(1,:) available(2,:)}, ...
                    {'parent',m2});
                % selections with first dim on y-axis, second dim on x-axis
                for i = 1:length(layout.x)
                    for j = 1:length(layout.y)
                        d = [layout.y(j) layout.x(i)];
                        available{1,i,j} = [header(d).dimID];
                        available{2,i,j} = fn_strcat({header(d).label},',');
                    end
                end
                fn_propcontrol(N,'selectiondimID', ...
                    {'menugroup' available(1,:) available(2,:)}, ...
                    {'parent',m2});
            end

            %             uimenu(m2,'label','ND selection...','separator','on', ...
            %                 'callback',@(u,e)set(N,'selectiondimID','prompt'))
            % (stop)
            if ~isempty(N.selectiondimID)
                uimenu(m,'label','Stop selection control', ...
                    'callback',@(u,e)set(N,'selectiondimID',[]))
            end

            % Options below should not be displayed if there is no active
            % selection control
            if isempty(N.selectiondimID), return, end

            % Selection options
            uimenu(m,'label','Clear selections','separator','on', ...
                'callback',@(u,e)N.selectionfilter.updateSelection('reset'))
            if length(N.selectiondimID) == 2
                fn_propcontrol(N,'selection2Dshape', ...
                    {'menuval' {'poly', 'free', 'rect', 'ellipse', 'ring', 'line', 'openpoly', 'freeline'}}, ...
                    'parent',m,'label','Shape');
            end
            if length(N.selectiondimID) == 1 && N.D.slice.header(seldim).ismeasure
                fn_propcontrol(N,'selectionround1Dmeasure','menu', ...
                    {'parent',m,'label','Round 1D selections to data indices','separator','on'});
            end
            fn_propcontrol(N,'selectionpromptname','menu', ...
                {'parent',m,'label','Prompt for name of new selections'});
            fn_propcontrol(N,'selectionshow', ...
                {'menuval', {'shape+name' 'shape' 'name' 'center'}}, ...
                {'parent',m,'label','Display mode','separator','on'});       
            % selection edit mode not ready yet!!! (need first to convert
            % selection from slice to graph)
            %             fn_propcontrol(N,'selectionedit','menu', ...
            %                 {'parent',m,'label','Selections modifyable'});

            % Load/save selections
            uimenu(m,'label','Load...','separator','on', ...
                'callback',@(u,e)N.selectionload())
            uimenu(m,'label','Save','enable',onoff(~isempty(N.selectionsavefile)), ...
                'callback',@(u,e)N.selectionsave(N.selectionsavefile))
            uimenu(m,'label','Save as...', ...
                'callback',@(u,e)N.selectionsave())
        end
        function selslice = selectionMouse(N, message)
            % function selslice = selectionMouse(N [,message])
            %---
            % 'pointaction' argument is a function handle that sets a
            % specific function to execute when drag action occurs to stay
            % on a point: this is used to raise a context menu when
            % right-clicking on a selection
            %
            % if second argument 'message' is set, this message will be
            % displayed under the mouse cursor, and it is assumed that the
            % first point of the selection was not clicked yet; otherwise
            % it is assumed that the first point of the selection was
            % already clicked
            
            seldim = N.selectiondim;
            selnd = length(seldim);

            % check dimension location
            dim_location = fn_strcat(N.D.layoutID.dim_locations(seldim),',');
            if ~ismember(dim_location, {'x' 'y' 'xy' 'x,y' 'y,x'})
                disp(['selection in location ''' dim_location ''' not handled'])
                selslice = [];
                return
            end
            
            % two different behaviors for the fn_mouse function
            if nargin<2
                % we assume the first point of the selection was already
                % clicked
                popt = '-';
                msgopt = {};
            else
                % we assume the first point of the selection was not
                % clicked yet, and a message will be displayed
                popt = '';
                msgopt = {message};
            end

            % define selection
            if selnd == 1
                % user interaction
                if isscalar(dim_location)
                    polyax = fn_mouse(N.ha,[dim_location 'segment' popt], msgopt{:});
                    switch dim_location
                        case 'x'
                            polyax = [polyax; 0 0];
                        case 'y'
                            polyax = [0 0; polyax];
                    end
                else
                    polyax = fn_mouse(N.ha,['rectaxp' popt], msgopt{:});
                end

                % convert to slice indices in the selected
                % dimension
                ijk0 = round(N.graph.graph2slice(polyax(:,1)));
                polyslice = N.graph.graph2slice(polyax, 'subdim', seldim, 'ijk0', ijk0);
                polyslice = sort(polyslice);

                % create selection in slice indices coordinates
                sz = N.D.slice.sz(seldim); % size of data in the dimension where selection is made
                if N.D.slice.header(seldim).categorical
                    selslice = xplr.selectionnd('indices',round(polyslice(1)):round(polyslice(2)),sz);
                elseif diff(polyslice)==0
                    selslice = xplr.selectionnd('point1D',polyslice(1));
                elseif N.selectionround1Dmeasure
                    selslice = xplr.selectionnd('line1D',round(polyslice)+[-.5 .5]);
                else
                    selslice = xplr.selectionnd('line1D',polyslice);
                end
            elseif selnd == 2
                % user interaction
                mouseselmode = fn_switch(N.selection2Dshape, ...
                    'line','segment',{'poly' 'openpoly'},'polypt','freeline','free', ...
                    'ellipse','ellipse*','ring','ring*', ...
                    N.selection2Dshape);
                seltype = fn_switch(N.selection2Dshape, ...
                    {'poly','free'},'poly2D','rect','rect2D', ...
                    'ellipse','ellipse2D','ring','ring2D', ...
                    {'line','openpoly','freeline'},'openpoly2D');
                polyax = fn_mouse(N.D.ha,[mouseselmode popt], msgopt{:});

                % create selection in graph coordinates
                selax = xplr.selectionnd(seltype,polyax);
                % if selection is too small, convert it to a single point
                selax.checkpoint(.005)
                
                % convert to slice coordinates
                selslice = N.graph.selection2slice(seldim,selax);
            end
        end
        function sel = get.selection(N)
            F = N.selectionfilter;
            if isempty(F)
                sel = [];
            else
                sel = F.selection;
            end
        end
        function set.selectiondimID(N,dimID)
            % function setselectiondimID(N,dimID)
            %---
            % Select the dimension(s) for which selections are displayed and
            % made in the display
            
            % special: prompt user for selecting dimensions
            if ischar(dimID) && strcmp(dimID,'prompt')
                header = N.D.slice.header;
                seldim = listdlg( ...
                    'PromptString', 'Select up to 2 dimensions', ...
                    'ListString', {header.label});
                dimID = [header(seldim(1:min(2,end))).dimID];
            end

            % check dimension(s)
            nd = length(dimID);
            if ~ismember(nd,[0 1 2])
                error 'number of dimension for selection display must be 1 or 2'
            end
            dim = N.D.slice.dimensionNumber(dimID);
            if iscell(dim), error 'some dimension is not present in the slice data', end
            singleton = (N.D.slice.sz(dim)==1);
            dimID(singleton) = []; % no selection in singleton dimension(s)
            dim(singleton) = [];
            if isequal(dimID, N.selectiondimID)
                return
            end
            
            % set property
            N.selectiondimID = dimID;
            
            % disconnect from previous filter
            F = N.selectionfilter;
            if ~isempty(F)
                xplr.bank.unregisterFilter(F,N)
                N.disconnect(F)
            end
            
            % no selection?
            if isempty(dimID)
                N.selectionfilter = [];
                N.displayselection()
                return
            end
            
            % find filter to connect to
            headerin = N.D.slice.header(dim); 
            F = xplr.bank.getFilterFilter(1, headerin, N); % xplr.filter object
            if isempty(F)
                % filter needs to be created
                error 'not implemented yet'
            end
            N.selectionfilter = F;
            
            % watch changes in filter to update display!
            N.addListener(F,'ChangedOperation',@(u,e)selectionfilterchange(N,e));
            
            % update display
            N.displayselection()
            
        end
        function checkselectionfilter(N)
            % Check whether current dimension for selections display is
            % still valid, i.e. whether the connected filter still fits a
            % dimension in the new slice
            if isempty(N.selectiondimID), return, end
            if ~all(ismember(N.selectiondimID,[N.D.slice.non_singleton_header().dimID]))
                N.selectiondimID = [];
            end
        end
        function selectionfilterchange(N,e)
            % Update selection display
            N.displayselection(e.flag,e.ind)
        end
        function displayselection(N,flag,ind)
            % function displayselection(N)
            % function displayselection(N,'reset')
            % function displayselection(N,'all')
            % function displayselection(N,'referentialchanged')
            % function displayselection(N,'new',ind)
            % function displayselection(N,'chg',ind)
            % function displayselection(N,'chg&new',{indchg indnew})
            % function displayselection(N,'chg&rm',{indchg indrm})
            % function displayselection(N,'remove',ind)
            % function displayselection(N,'perm',perm)
           
            if isempty(N.selectiondimID)
                deleteValid(N.selectiondisplay)
                N.selectiondisplay = [];
                return
            end
            
            % before all, check that selections can be displayed
            selectionaxis = [N.D.layoutID.dim_locations{N.selectiondim}];
            if ~ismember(selectionaxis, {'x' 'y' 'xy' 'yx'})
                disp('selections cannot be displayed')
                deleteValid(N.selectiondisplay)
                N.selectiondisplay = [];
                return
            end
            
            % or display update
            if nargin<2, flag = 'all'; end
            if nargin<3, ind = 1:length(N.selection); end
            switch flag
                case {'all' 'reset' 'new'}
                    % delete current display
                    if fn_ismemberstr(flag,{'all','reset'})
                        deleteValid(N.selectiondisplay)
                        N.selectiondisplay = [];
                    end
                    % draw new selections
                    for idx=ind, displayonesel(N,idx,'new'); end
                    % keep cross above selections, with the cross center
                    % N.cross(3) at the very top
                    try uistack(N.cross([3 1 2]),'top'), end % can fail when D.resetDisplay is invoked
                case {'add','chg', 'affinity'}
                    % might be several indices
                    for k=ind, displayonesel(N,k,'pos'); end
                case 'referentialchanged'
                    % it is not the positions of selections that have
                    % changed, but the referential of these positions
                    % relative to the main display axes: we need then to
                    % recompute the positions of each selection display
                    % inside this new referential
                    for k = 1:length(N.selection), displayonesel(N, k, 'pos'), end
                case 'remove'
                    % delete selected selections
                    deleteValid(N.selectiondisplay(ind))
                    N.selectiondisplay(ind) = [];
                    % index (and therefore displayed name) of some other selections have changed
                    for k = min(ind):length(N.selection), displayonesel(N, k, 'name'), end
                case 'perm'
                    perm = ind;
                    N.selectiondisplay = N.selectiondisplay(perm);
                    % index (and therefore displayed name) of some selections have changed
                    for k = find(perm~=1:length(N.selection)), displayonesel(N, k, 'name'), end
                otherwise
                    error('unknown flag ''%s''', flag)
            end
                
        end
        function displayonesel(N,k,flag,varargin)
            % function displayonesel(D,k,'new')       - selection at index k is new
            % function displayonesel(D,k,'pos')       - position has changed
            % function displayonesel(D,k,'name')      - name or color has changed
            
            % flags: what kind of update
            [flagnew, flagpos, flagname] = fn_flags('new','pos','name',flag);
            
            % selection in slice coordinates
            selij = N.selection(k);
            name = N.selectionfilter.headerout.getItemNames{k};
            
            % convert selection to displayable polygon
            seldim = N.selectiondim;
            if flagnew || flagpos
                selectionaxis = [N.D.layoutID.dim_locations{seldim}];
                if ~ismember(selectionaxis, {'x' 'y' 'xy' 'yx'})
                    error('selection cannot be displayed')
                end
                % selection shape
                [polygon, center] = N.D.graph.selectionMark(seldim,selij);
                % additional points for when in selection edit mode
                if N.selectionedit
                    switch selij.type
                        case {'poly2D','mixed','point2D','line2D'} % TODO: not sure about 'point2D'
                            shape_edit_poly = polygon(:,1:end-1); % the last point is a repetition of the 1st one
                            shape_edit_info = [];
                        case 'rect2D'
                            shape_edit_poly = polygon(:,1:4); % the 5th point of polygon is a repetition of the 1st one
                            shape_edit_info = [selij.shapes.points' selij.shapes.vectors'];
                        case {'ellipse2D' 'ring2D'}
                            c = selij.shapes.points;
                            u = selij.shapes.vectors;
                            e = selij.shapes.logic;
                            shape_edit_poly = [c-u c+u];
                            shape_edit_info = {c u e};
                        otherwise
                            error programming
                    end
                end
            end
            
            % color
            if flagnew || flagname
                colors = fn_colorset;
                col = colors(mod(k-1,size(colors,1))+1,:);
                namerotation = fn_switch(isscalar(N.selectiondimID) && length(name)>3, 45, 0);
            end
            
            % Create / update objects
            if flagnew
                hl = struct('shape',[],'label',[],'cross',[],'handles',[]);
                % selection shape
                if strfind(N.selectionshow,'shape')
                    hl.shape = line(polygon(1,:),polygon(2,:),'Parent',N.D.ha);
                end
                % name
                if strfind(N.selectionshow,'name')
                    hl.label = text(center(1),center(2),name,'Parent',N.D.ha, ...
                        'horizontalalignment','center','verticalalignment','middle', ...
                        'rotation',namerotation);
                end
                % center marked with a cross
                if strfind(N.selectionshow,'center')
                    hl.cross = line(center(1),center(2),'Parent',N.D.ha, ...
                        'LineStyle','none','Marker','+','MarkerSize',4);
                end
                % handles to modify selection
                if N.selectionedit
                    hl.handles = line(shape_edit_poly(1,:),shape_edit_poly(2,:),'Parent',N.D.ha, ...
                        'LineStyle','none','marker','.','UserData',shape_edit_info);
                end
                set(struct2array(hl),'Color',col, ...
                    'ButtonDownFcn',@(u,e)N.selection_clicked(u));
                if isempty(N.selectiondisplay)
                    if k~=1, error 'attempting to create N.selectiondisplay at initial index different from 1', end
                    N.selectiondisplay = hl;
                else
                    N.selectiondisplay(k) = hl;
                end            
            else
                hl = N.selectiondisplay(k);
                if flagpos
                    set(hl.label,'position',center)
                    set(hl.shape,'xdata',polygon(1,:),'ydata',polygon(2,:))
                    set(hl.cross,'xdata',center(1),'ydata',center(2))
                    if N.selectionedit
                        set(hl.handles,'xdata',shape_edit_poly(1,:),'ydata',shape_edit_poly(2,:), ...
                            'UserData',shape_edit_info);
                    end
                elseif flagname
                    set(hl.label,'string',name,'rotation',namerotation)
                    set(struct2array(hl),'color',col)
                end
            end
        end
        function set.selectionshow(N,value)
            N.selectionshow = value;
            % we must see the shape when selection edit is active
            if N.selectionedit && ~any(strfind(N.selectionshow,'shape'))
                N.selectionedit = false;
            else
                N.displayselection('all')
            end
        end
        function set.selectionedit(N,value)
            N.selectionedit = value;
            % we must see the shape when selection edit is active
            if N.selectionedit && ~any(strfind(N.selectionshow,'shape'))
                N.selectionshow = 'shape+name';
            else
                N.displayselection('all')
            end
        end
        function selection_clicked(N,u)
            % click on selection u: 
            % - if selection edit mode is active, some specific actions can
            %   happen
            % - if right-click and no mouse motion, raise a context menu
            % - if none of the two types of actions above are triggered,
            %   execute normal axes callback
            
            % first get index of the clicked selection
            idx = fn_find(@(hl)any(struct2array(hl)==u),N.selectiondisplay,'first');
            
            click_type = get(N.D.V.hf,'SelectionType');
            if strcmp(click_type,'alt')
                % perform a regular new selection, but if mouse did not
                % move, show the context menu instead of creating a point
                % selection
                selslice = N.selectionMouse();
                if ~ispoint(selslice)
                    % prompt for selection name
                    options = {};
                    if N.selectionpromptname
                        name = inputdlg('Selection name','xplor');
                        if ~isempty(name), options = {'Name' name{1}}; end
                    end

                    % update filter
                    N.selectionfilter.updateSelection('new',selslice,options{:})
                else
                    % instead of creating a point selection, raise the
                    % context menu
                    N.selection_context_menu(idx)
                end
            elseif N.selectionedit
                flag = '';
                if u == N.selectiondisplay{idx}(1)
                    flag = 'line';
                elseif u == N.selectiondisplay{idx}(end)
                    flag = 'handle';
                end
                if ~isempty(flag)
                    N.selection_edit(N,idx,flag)
                else
                    N.axes_click()
                end
            else
                N.axes_click()
            end
        end
        function selection_edit(N,idx,flag)
            % function selection_edit(N,idx,'line|handle|redraw')
            
            switch flag
                case 'redraw'
                    selslice = N.selectionMouse('select new shape');
                    N.selectionfilter.updateSelection('chg',idx,selslice)
            end
            
            % code below is not ready
            return
            
            htl = N.selectiondisplay{idx};
            hshape = hl(1);
            hhandle = hl(end);
            hl = htl(2:end);
            [flagpt, flaglin] = fn_flags({'handle','line'},flag);
            p = get(N.ha,'currentpoint'); p = p(1,1:2)';
            polymark = [get(hl(2),'xdata'); get(hl(2),'ydata')];
            seldimsnum = N.seldims-'w';
            selectionmarks = N.SI.selection.getselset(seldimsnum).singleset;
            switch get(N.hf,'selectiontype')
                case 'normal'               % MOVE POINT
                    shapetype = selectionmarks(idx).type;
                    switch shapetype
                        case {'poly2D' 'mixed' 'line2D' 'openpoly2D'}
                            % note that first and last point in polygon are
                            % the same!
                            if flagpt
                                % closest point
                                dist = sum(fn_add(polymark,-p).^2);
                                [dum idx] = min(dist); idx = idx(1); %#ok<*ASGLU>
                                if idx==1 && ~fn_ismemberstr(shapetype,{'line2D' 'openpoly2D'})
                                    % need to move both the first and last point (which is a repetition of the first)
                                    idx=[1 size(polymark,2)]; 
                                end 
                            else
                                % closest segment (in fact, closest line)
                                a = polymark(:,1:end-1);
                                b = polymark(:,2:end);
                                ab = b-a;
                                ab2 = sum(ab.^2);
                                ap = fn_add(p,-a);
                                abap = ab(1,:).*ap(2,:)-ab(2,:).*ap(1,:);
                                dist = abs(abap) ./ ab2;
                                [dum idx] = min(dist); idx = idx(1);
                                polymark = [a(:,1:idx) p b(:,idx:end)];
                                set(hl,'xdata',polymark(1,:),'ydata',polymark(2,:))
                                idx = idx+1;
                            end
                            fn_moveobject(hl,'point',idx)
                            seleditupdateslice(N,idx)
                        case 'rect2D'
                            desc = getappdata(hl(2),'description');
                            x = desc(1); y = desc(2);
                            w = desc(3); h = desc(4);
                            x2 = x+w; y2 = y+h;
                            if flagpt
                                % move corner
                                pol = [x x2 x2 x; y y y2 y2]; % anti-clockwise from (x,y)
                                dist = sum(fn_add(pol,-p).^2);
                                [dum idx] = min(dist); idx = idx(1);
                            else
                                % move edge
                                dist = abs([p(2)-y p(1)-x2 p(2)-y2 p(1)-x]);
                                [dum idx] = min(dist);
                            end
                            col = get(hl(1),'color');
                            set(hl,'color',.5*[1 1 1])
                            chgrectangle(N.ha,hl,flagpt,idx,desc)
                            fn_buttonmotion({@chgrectangle,N.ha,hl,flagpt,idx,desc},N.hf);
                            set(hl,'color',col)
                            seleditupdateslice(N,idx)
                        case {'ellipse2D' 'ring2D'}
                            desc = getappdata(hl(2),'description');
                            if flagpt
                                % closest of two anchor points
                                dist = sum(fn_add(polymark,-p).^2);
                                [dum idx] = min(dist); idx = idx(1);
                            elseif strcmp(shapetype,'ellipse2D')
                                % eccentricity
                                idx = 0;
                            else
                                % eccentricity or secondary radius?
                                polygon = [get(hl(1),'xdata'); get(hl(1),'ydata')];
                                dist = sum(fn_add(polygon,-p).^2);
                                [dum idx] = min(dist);
                                if idx<length(polygon)/2
                                    % eccentricity
                                    idx = 0;
                                else
                                    % secondary radius
                                    idx = -1;
                                end
                            end
                            col = get(hl(1),'color');
                            set(hl,'color',.5*[1 1 1])
                            chgellipse(N.ha,hl,idx,desc)
                            fn_buttonmotion({@chgellipse,N.ha,hl,idx,desc},N.hf);
                            set(hl,'color',col)
                            seleditupdateslice(N,idx)
                        otherwise
                            error programming
                    end
                case 'extend'               % MOVE SHAPE
                    if flagpt
                        dp = fn_moveobject(htl);
                        seleditupdateslice(N,idx,dp)
                    elseif flaglin
                        dp = fn_moveobject([N.seldisp{:}]);
                        seleditupdateslice(N,1:length(N.seldisp),dp) % move all shapes
                    end
                case 'alt'                  % REMOVE
                    if fn_ismemberstr(selectionmarks(idx).type, ...
                            {'poly2D','mixed'}) && flagpt
                        % closest point -> remove vertex
                        dist = sum(fn_add(polymark,-p).^2);
                        [dum idx] = min(dist); idx = idx(1);
                        if idx==1, idx=[1 size(polymark,2)]; end
                        polymark(:,idx) = [];
                        selax = selectionND('poly2D',polymark);
                        sel = AX2IJ(N.SI,selax);
                        updateselection(N.SI,'change',idx,sel)
                    else
                        % replace the whole shape
                        set(hl,'visible','off')
                        mouseselmode = fn_switch(N.shapemode, ...
                            {'poly' 'openpoly'},'polypt','freeline','free', ...
                            N.shapemode);
                        TYPE = fn_switch(N.shapemode,{'poly','free'},'poly2D','rect','rect2D', ...
                            'ellipse','ellipse2D','ring','ring2D', ...
                            'segment','line2D',{'openpoly','freeline'},'openpoly2D');
                        polyax = fn_mouse(N.ha,mouseselmode,'select new shape');
                        selax = selectionND(TYPE,polyax);
                        sel = AX2IJ(N.SI,selax);
                        updateselection(N.SI,'change',idx,sel)
                        set(hl,'visible','on')
                    end
            end
        end
        function selection_context_menu(N,idx)
            % init context menu
            if isempty(N.selection_context)
                % create context menu for the first time
                N.selection_context = uicontextmenu(N.D.V.hf);
            end
            m = N.selection_context;
            delete(get(m,'children'))
            
            % remove selection
            uimenu(m,'label','remove selection', ...
                'callback',@(u,e)N.selectionfilter.updateSelection('remove',idx))
            
            % replace selection
            uimenu(m,'label','re-draw selection', ...
                'callback',@(u,e)N.selection_edit(idx,'redraw'))
            
            % change selection number
            nsel = N.selectionfilter.nsel;
            if idx ~= 1
                uimenu(m,'label','make this selection first', ...
                    'callback',@(u,e)N.selectionfilter.updateSelection('perm',[idx setdiff(1:nsel,idx)]))
            end
            if idx ~= nsel
                uimenu(m,'label','make this selection last', ...
                    'callback',@(u,e)N.selectionfilter.updateSelection('perm',[setdiff(1:nsel,idx) idx]))
            end
            
            % make menu visible
            p = get(N.D.V.hf,'currentpoint'); p = p(1,1:2);
            set(m,'Position',p,'Visible','on')
        end
        function selectionsave(N,fname)
            if nargin<2 || isempty(fname)
                prompt = 'Select file for saving selections';
                fname = fn_savefile('*.xpls',prompt,N.selectionsavefile);
                if isequal(fname,0), return, end
            end
            N.selectionfilter.savetofile(fname);
            N.selectionsavefile = fname;
        end
        function selectionload(N,fname)
            if nargin<2
                fname = fn_getfile('*.xpls','Select selections file'); 
                if isequal(fname,0), return, end
            end
            try
                N.selectionfilter.loadfromfile(fname);
                N.selectionsavefile = fname;
            catch ME
                errordlg(ME.message)
            end
        end
   end

    % Slider and scroll wheel callbacks: change zoom
    methods
        function chgzoom(N,f,obj)
            dim = N.D.activedim.(f);
            if isempty(dim), return, end
            % linked object
            Z = N.D.zoomfilters(dim);
            % prevent unnecessary update of slider display
            c = disableListener(N.sliders.(f));
            % set value
            if isequal(obj.value,obj.minmax)
                setZoom(Z,':')
            else
                setZoom(Z,obj.value)
            end
        end
        function Scroll(N,nscroll)
            p = get(N.D.ha,'currentpoint'); p = p(1,1:2);
            origin = row(N.graph.graph2slice(p)); % current point in data coordinates
            zoomfactor = 1.5^nscroll;
            dim = [N.D.activedim.x N.D.activedim.y];
            if isempty(dim), return, end
            % This commented code had been put to replace the line below,
            % but it seems that the effect is less intuitive. Let's go back
            % to the previous code and see if we get errors or unintuitive
            % behaviors to decide what to do. (TD 12/11/2019)
            %             if nscroll<0 && ~any([N.D.layout.xy N.D.layout.yx])
            %                 % it does not make sense to zoom-in in a dimensions which
            %                 % does not fill its available space due to aspect ratio
            %                 % constraints
            %                 dim(N.graph.filling(dim)<1) = [];
            %             end
            %             zoom = N.graph.getZoom(dim); %,'effective');
            zoom = N.graph.getZoom(dim,'effective');
            newzoom = fn_add(origin(dim), fn_mult(zoomfactor,fn_subtract(zoom,origin(dim))));
            %fprintf('%.2f -> %.2f\n',diff(zoom),diff(newzoom))
            N.D.zoomslicer.setZoom(dim,newzoom)
        end
    end

    % Update upon changes in active dim and zoom
    methods
        function connectZoomFilter(N,f)
            % both x and y?
            if nargin<2
                connectZoomFilter(N,'x')
                connectZoomFilter(N,'y')
                return
            end
            % slider object and corresponding data dimension
            obj = N.sliders.(f);
            d = N.D.activedim.(f);
            % disconnect from previous zoomfilters
            Zold = N.zoomfilters.(f);
            if ~isempty(Zold), N.disconnect(Zold), end
            % no active dim?
            if isempty(d)
                set(obj,'visible','off')
                return
            end
            % update slider display to reflect zoom in the specified
            % dimension
            Z = N.D.zoomfilters(d);
            N.zoomfilters.(f) = Z;
            set(obj,'visible','on','minmax',[.5 Z.headerin.n+.5],'value',Z.zoomvalue)
            % watch changes in zoom
            function zoomchange(u,e)
                if strcmp(e.type,'zoom')
                    set(obj,'value',Z.zoomvalue)
                end
            end
            N.addListener(Z,'ChangedOperation',@zoomchange);
        end
    end
    
    % Tool used by both cross and selection display: is point out of display
    methods
        % return true if point (graph coordinates) is part of the slice
        % data displayed by converting the point to slice coordinates and
        % test if its between minimal and maximal values for all dimensions
        %
        % @param point: 2xn double 
        % @param perdim: whether to return a single logical value per point
        %                (point is/isn't inside slice, default behavior)
        %                or a vector (for each dimension, if perdim=true)
        % @return output: 1xn boolean
        function output = isPointOutOfDisplay(N, point, perdim)
            if nargin<3, perdim = false; end
            
            % get slice indices corresponding to the point
            ijk = N.graph.graph2slice(point,'invertible',true)'; % n x ndim
            
            output = isIndexOutOfDisplay(N, ijk, perdim);
        end
        function output = isIndexOutOfDisplay(N, ijk, perdim) 
            % input
            if nargin<3, perdim = false; end
            if isvector(ijk) && size(ijk,2) ~= N.D.slice.nd
                ijk = row(ijk);
            end
            
            % get the min and max slice values of the data displayed
            zoomSliceValues = N.graph.getZoom('displaylimit'); % 2 x ndim
            
            % set the ouput to zeros (they will be set to one if one of the
            % dimension if it's out of display)
            if perdim
                output = bsxfun(@lt,ijk,zoomSliceValues(1,:)) ...
                    | bsxfun(@gt,ijk,zoomSliceValues(2,:)); % n x ndim
            else
                output = min(ijk,1)<zoomSliceValues(1,:) ...
                    | max(ijk,1)>zoomSliceValues(2,:); % 1 x ndim
            end
        end
    end
    
    
    
end