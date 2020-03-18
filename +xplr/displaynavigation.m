classdef displaynavigation < xplr.graphnode
% display navigation

    properties (SetAccess='private')
        D                                       % parent xplr.viewdisplay
        ha
        hf
        graph
        crossCenter
        cross                                   % display cross selector
        sliders = struct('x',[],'y',[]);        % slider objects
        zoomfilters = struct('x',[],'y',[]);    % connected zoom filters
        dimfilters = {};
    end
    properties
        selectionfilter         % filter being modified by the selections
        selectiondisplay        % displays of selectionnd
    end
    properties (SetObservable)
        selectiondimID          % dimensions to which these selections apply, identified by its id
        selection2Dshape = 'ellipse'; % 'xsegment', 'ysegment', 'poly', 'free', 'rect', 'ellipse', 'ring', 'segment', 'openpoly', 'freeline'
        selectionround1Dmeasure = true; 
        selectionadvanced = false
    end
    properties (Dependent, SetAccess='private')
        layout
        selection               % list of selectionnd object from the filter
        selectiondim            % dimensions to which these selections apply, identified by its number
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
            
            % connect sliders to the active dimensions of the display
            % (note that this is in fact redundant with call in
            % viewdisplay.slicechange when viewdisplay object is created)
            connectZoomFilter(N)

            % mouse actions
            set(D.ha,'buttondownfcn',@(u,e)Mouse(N))

            % scroll wheel zooming
            fn_scrollwheelregister(D.ha,@(n)N.Scroll(n))
            
            % selection menu
            uimenu(N.hf,'Label','Selection','callback',@(m,e)N.selectionMenu(m))
        end
        function init_buttons(N)
        % function init_buttons(N)
        % 3 buttons that control clipping
            
            % first button to adjust clipping with mouse movements:
            % display image on it indicating how image luminance and
            % contrast change upon mouse movements
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
    end
    
    % Get dependent
    methods
        function layout = get.layout(N)
            layout = N.D.layout;
        end
        function d = get.selectiondim(N)
            d = N.D.slice.dimensionByID(N.selectiondimID);
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
                        case 'mean(line)'
                            if strcmp(N.D.displaymode,'time courses'), clipcenter = 0; end
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
    
    % Mouse actions
    methods
        function Mouse(N, flag)
            pointonly = (nargin==2 && strcmp(flag,'pointonly'));
            point =  get(N.D.ha,'CurrentPoint'); point = point(1,[1 2])';
            activedim = [N.D.activedim.x N.D.activedim.y];
            switch get(N.hf,'SelectionType')
                case 'normal'
                    % zoom in or select point
                    if pointonly                        
                        dozoom = false;
                    else
                        % rect = fn_mouse(N.ha,'rectangle-');
                        rect = fn_mouse(N.ha,'rectangle-');
                       % rect=[-0.3911,-0.3911,-0.2189,-0.2189;0.2391,0.0468,0.0468,0.2391];
                        dozoom = any(any(abs(diff(rect,1,2))>1e-2));
                    end
                    if dozoom
                        ijk = N.graph.graph2slice(rect(:,[1 3]));
                        zoom = ijk(activedim,:)';
                        for i=1:length(activedim), zoom(:,i) = sort(zoom(:,i)); end
                        N.D.zoomslicer.setZoom(activedim,zoom)
                    else
                        N.manualclickmovecross(point);
                    end
                case 'open'
                    % zoom reset
                    zoom = repmat(':',1,length(activedim));
                    N.D.zoomslicer.setZoom(activedim,zoom)
                case 'alt'
                    % on right click: create a new selection in N.selection
                    % depending on the parameter selectionshape
                    
                    % TODO: the code below is work in progress, very
                    % specific. It has to be generalized.
                    
                    d = N.selectiondim;
                    switch length(d)
                        
                        case 1
                            
                            dim_location = N.layout.dim_locations{d};
                            if ~ismember(dim_location, {'x' 'y' 'xy'})
                                disp(['selection in location ''' dim_location ''' not handled'])
                                return
                            end
                            
                            % user interaction
                            if isscalar(dim_location)
                                segment = fn_mouse(N.ha,[dim_location 'segment-']);
                                switch dim_location
                                    case 'x'
                                        segment = [segment; 0 0];
                                    case 'y'
                                        segment = [0 0; segment];
                                end
                            else
                                segment = fn_mouse(N.ha,'rectaxp-');
                            end
                            
                            % convert to slice indices in the selected
                            % dimension
                            ijk = N.graph.graph2slice(segment);
                            isegment = sort(ijk(d,:));
                            
                            % create selection in slice indices coordinates
                            sz = N.D.slice.sz(d); % size of data in the dimension where selection is made
                            if N.D.slice.header(d).categorical
                                oneDselection = xplr.selectionnd('indices',round(isegment(1)):round(isegment(2)),sz);
                            elseif diff(isegment)==0
                                oneDselection = xplr.selectionnd('point1D',isegment(1));
                            elseif N.selectionround1Dmeasure
                                oneDselection = xplr.selectionnd('line1D',round(isegment)+[-.5 .5]);
                            else
                                oneDselection = xplr.selectionnd('line1D',isegment);
                            end
                            
%                             % update selection
%                             if isempty(N.selection)
%                                N.selection=oneDselection; 
%                             else
%                                N.selection(end+1)=oneDselection;
%                             end
                                                        
                            % update filter
                            N.selectionfilter.updateSelection('new',oneDselection)
                            
                            % update display
                            N.displayselection('new',length(N.selection))
                        case 'poly'
                            poly = fn_mouse(N.ha,'poly');
                            ijk = N.graph.graph2slice(poly,...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            
                            polySelection = xplr.selectionnd('poly2D',ijk(d,:));
                             if isempty(N.selection)
                               N.selection=polySelection; 
                            else
                                N.selection(end+1)=polySelection;
                            end
                            N.displayselection('new',length(N.selection))
                        case 'free'
                            poly = fn_mouse(N.ha,'free');
                            ijk = N.graph.graph2slice(poly,...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            
                            polySelection = xplr.selectionnd('poly2D',ijk(d,:));
                             if isempty(N.selection)
                               N.selection=polySelection; 
                            else
                                N.selection(end+1)=polySelection;
                            end
                            N.displayselection('new',length(N.selection))
                        case 'rect'
                            rectangle = fn_mouse(N.ha,'rectangle-');
                            
                            ijk1 = N.graph.graph2slice(rectangle(:,1),...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            ijk2 = N.graph.graph2slice(rectangle(:,3),...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            vector=ijk1-ijk2;
                            width = abs(vector(d(1),1));
                            height= abs(vector(d(2),1));
                            x = min(ijk1(d(1),1),ijk2(d(1),1));
                            y = min(ijk1(d(2),1),ijk2(d(2),1));
                            rectangleSelection = xplr.selectionnd('rect2D',[x,y,width,height]);
                            if isempty(N.selection)
                               N.selection=rectangleSelection; 
                            else
                                N.selection(end+1)=rectangleSelection;
                            end
                            N.displayselection('new',length(N.selection))
                        case 'ellipse'
                            ellipse = fn_mouse(N.ha,'ellipse-*');
        %                     disp("ellipse");
        %                     ellipse{2} = ellipse{1} + ellipse{2};
        %                     ellipse{1}
        %                     ellipse{2}
        %                     ellipse{3}
        %                     disp("end ellipse");
                            ijk1 = N.graph.graph2slice(ellipse{1},...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            ijk2 = N.graph.graph2slice(ellipse{2},...
                                'mode','vector', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));

                            ijk = {ijk1(d),ijk2(d), ellipse{3}};
                            ellipseSelection = xplr.selectionnd('ellipse2D',ijk);

                            if isempty(N.selection)
                               N.selection=ellipseSelection; 
                            else
                                N.selection(end+1)=ellipseSelection;
                            end

                            N.displayselection('new',length(N.selection))
                        case 'ring'
                            ring = fn_mouse(N.ha,'ring*');
        %                     disp("ellipse");
        %                     ellipse{2} = ellipse{1} + ellipse{2};
        %                     ellipse{1}
        %                     ellipse{2}
        %                     ellipse{3}
        %                     disp("end ellipse");
                            ijk1 = N.graph.graph2slice(ring{1},...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            ijk2 = N.graph.graph2slice(ring{2},...
                                'mode','vector', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));

                            ijk = {ijk1(d),ijk2(d), [ring{3},ring{4}]};
                            ellipseSelection = xplr.selectionnd('ring2D',ijk);

                            if isempty(N.selection)
                               N.selection=ellipseSelection; 
                            else
                                N.selection(end+1)=ellipseSelection;
                            end

                            N.displayselection('new',length(N.selection))
                        case 'segment'
                            poly = fn_mouse(N.ha,'segment');
                            ijk = N.graph.graph2slice(poly,...
                                'mode','point', ...
                                'subdim',d', ...
                                'ijk0',ones(N.D.nd,1));
                            
                            polySelection = xplr.selectionnd('poly2D',ijk(d,:));
                             if isempty(N.selection)
                               N.selection=polySelection; 
                            else
                                N.selection(end+1)=polySelection;
                            end
                            N.displayselection('new',length(N.selection))
                        otherwise
                            % don't exist yet
                        
                    end
            end
        end
    end
    
    % Cross point selection
    methods
        function connectPointFilter(N,dim,key)
            if nargin < 3
               key = 1; 
            end

            if nargin < 2
                disconnectPointFilter(N)
                dim = 1:N.D.slice.nd;
            else
                disconnectPointFilter(N,dim)
            end
            
            for d = dim
                linkkey = key;
                head = N.D.slice.header(d);
                % no interest in creating and controling a filter for a
                % dimension with only 1 value
                if head.n ==1 
                    N.dimfilters{d} = [];
                    continue
                end
                % get filter from bank or create one for the header in this
                % dimension
                doshow=false;
                F = xplr.bank.getFilter(linkkey,head,doshow,N); % FilterAndPoint filter
                N.dimfilters{d} = F;
                % listen to the point filter event
                P = F.P; % Point filter
                N.addListener(P,'ChangedPoint',@(u,e)movedPoint(N))
            end
        end
        function disconnectPointFilter(N,dim)
            if nargin < 2
                dim = 1:length(N.dimfilters);
            end
            for d = dim
                F = N.dimfilters{d};
                if isempty(F), continue, end
                xplr.bank.unregisterFilter(F,N)
                N.disconnect(F.P)  % this is the same as F.disconnect(N)!
            end
        end
        function movedPoint(N)
            ijk = getPointIndexPosition(N);
            N.crossCenter = N.graph.slice2graph(ijk);
            update_cross_visibility(N);
        end
        function ijk = getPointIndexPosition(N)
            nd = N.D.slice.nd;
            ijk = ones(nd, 1);
            for d = 1:nd
                F = N.dimfilters{d};
                if isempty(F), continue, end
                ijk(d) = F.P.index0;
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
            N.cross(3) = line('Parent',N.D.ha,'xdata',0,'ydata',0,'marker','.','linestyle','none'); % a single point
            set(N.cross,'Color','k')
            
            %fn4D_dbstack
            %ij2 = D.SI.ij2;
            % scaling and translation
            %pt = IJ2AX(D.SI,ij2);
            N.crossCenter = [0 0];
            
            for i=1:3
                set(N.cross(i),'buttondownfcn',@(u,e)manualmovecross(N,i))
            end
        end
        function set.crossCenter(N, crossCenter)
            % set the property
          
            N.crossCenter = crossCenter;

            % move the cross
            set(N.cross(1),'XData',crossCenter([1 1]))
            set(N.cross(2),'YData',crossCenter([2 2]))
            set(N.cross(3),'XData',crossCenter(1),'YData',crossCenter(2))

        end
        function manualmovecross(N,il)
            if ~strcmp(get(N.hf,'selectiontype'),'normal')
                % not a left click: execute callback for axes
                Mouse(N)
                return
            end
            set(N.hf,'pointer',fn_switch(il,1,'left',2,'top',3,'cross'))
            anymove = fn_buttonmotion(@movecrosssub,N.hf,'moved?');
            set(N.hf,'pointer','arrow')
            if ~anymove
                % execute callback for axes
                Mouse(N, 'pointonly')
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
                manualclickmovecross(N,point)
            end
            
        end
        
        function manualclickmovecross(N,point)
            % move the cross to the selected point

            ijk = N.graph.graph2slice(point,'invertible',true);
            % update the point filters (only for dimensions where the point
            % remains within the slice)
            for d = find(~isOutOfDisplay(N,point,true))
                F = N.dimfilters{d};
                if ~isempty(F)
                    F.P.index = ijk(d);
                end
            end
        end
        
        function removecross(N)
            delete(N.cross)
        end      
        function update_cross_visibility(N)
           
            % if the slice has only one value in the dimension
            % displayed in abscisse then hide vertical bar
            x_singleton = true;
            x_isOutOfDisplay = false;
            ijk = getPointIndexPosition(N);
            zoom = N.graph.getZoom('displaylimit');
            
            for dimension = N.D.layout.x
                dimension_isSingleton = (N.D.V.slice.header(dimension).n == 1);
                if ~dimension_isSingleton, x_singleton = false; end
                
                % for all dimensions in layout.x, check if the crossCenter is
                % out of the display of this dimension. If the crossCenter
                % is out of display, the vertical bar will be hidden
                if ijk(dimension)<zoom(1,dimension) || ijk(dimension)>zoom(2,dimension)
                    x_isOutOfDisplay = true;
                end
            end
            if ~isempty(N.layout.griddisplay) && N.D.V.slice.header(N.layout.griddisplay).n > 1
                x_singleton = false;
            end
            
            % Hide the vertical if all dimensions on x are singletons or if
            % crossCenter is out of display on one dimension on x
            set(N.cross(1),'Visible',fn_switch(~(x_singleton|x_isOutOfDisplay)));
            
            % same things for horizontal bar
            y_singleton = true;
            y_isOutOfDisplay = false;
            for dimension = N.D.layout.y
                dimension_isSingleton = (N.D.V.slice.header(dimension).n == 1);
                if ~dimension_isSingleton, y_singleton = false; end

                if ijk(dimension)<zoom(1,dimension) || ijk(dimension)>zoom(2,dimension)
                    y_isOutOfDisplay = true;
                end
            end
            if ~isempty(N.layout.griddisplay) && N.D.V.slice.header(N.layout.griddisplay).n > 1
                y_singleton = false;
            end

            set(N.cross(2),'Visible',fn_switch(~(y_singleton|y_isOutOfDisplay)));
            
            updateCrossCenterVisibility(N);
        end
    end
    
    % ROI selection
    methods
        function selectionMenu(N,m)
            % Selection menu is populated when being activated
            delete(get(m,'children'))
            
            slice = N.D.slice;
            fn_propcontrol(N,'selectiondimID', ...
                {'menuval' {slice.header.dimID []} {slice.header.label}}, ...
                'parent',m,'label','in dimension(s)');
            if length(N.selectiondimID) == 2
                fn_propcontrol(N,'selection2Dshape', ...
                    {'menuval' {'ysegment', 'xsegment','poly', 'free', 'rect', 'ellipse', 'ring', 'segment', 'openpoly', 'freeline'}}, ...
                    'parent',m,'label','shape');
            end
            if ~isempty(N.selectiondimID)
                uimenu(m,'label','clear selections','separator','on', ...
                    'callback',@(u,e)N.selectionfilter.updateSelection('reset'))
            end
            if length(N.selectiondimID) == 1 && N.D.slice.header(N.selectiondim).ismeasure
                fn_propcontrol(N,'selectionround1Dmeasure','menu', ...
                    {'parent',m,'label','round selections to data indices','separator','on'});
                nextsep = 'off';
            else
                nextsep = 'on';
            end
            fn_propcontrol(N,'selectionadvanced','menu', ...
                {'parent',m,'label','advanced selection','separator',nextsep});
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
            
            % check dimension(s)
            nd = length(dimID);
            if ~ismember(nd,[0 1 2])
                error 'number of dimension for selection display must be 1 or 2'
            end
            dim = N.D.slice.dimensionByID(dimID);
            if iscell(dim), error 'some dimension is not present in the slice data', end
            singleton = (N.D.slice.sz(dim)==1);
            dimID(singleton) = []; % no selection in singleton dimension(s)
            dim(singleton) = [];
            if isequal(dimID, N.selectiondimID)
                return
            end
            N.selectiondimID = dimID;
            
            % disconnect from previous filter
            F = N.selectionfilter;
            if ~isempty(F)
                N.disconnect(F)
            end
            
            % no selection?
            if isempty(dimID)
                N.selectionfilter = [];
                return
            end
            
            % find filter to connect to
            headerin = N.D.slice.header(dim); 
            F = xplr.bank.getFilter(1, headerin); % xplr.filterAndPoint object
            if isempty(F)
                % filter needs to be created
                error 'not implemented yet'
            end
            F = F.F; % xplr.filter object
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
            if isempty(N.D.slice.dimensionByID(N.selectiondimID))
                N.selectiondimID = [];
            end
        end
        function selectionfilterchange(N,e)
            % Update selection display
            N.displayselection()
        end
        function displayselection(N,flag,ind,value)
            % @param flag: string 'all', 'new'
            % @param ind: integer
            % @param value:
            %
            % @return:
           
            if isempty(N.selectiondimID)
                deleteValid(N.selectiondisplay)
                N.selectiondisplay = [];
                return
            end
            
%             disp(['selection in dimensions ' num2str(N.selectiondim) ':'])
%             disp(N.selection)
%             disp(' ')
       

%             %N.selection(1).;
%             if ~isempty(N.selectiondisplay)
%                 deleteValid(N.selectiondisplay{:})
% 
%             end
%             
%             for i = 1:length(N.selection)
%                 N.displayonesel(i,'new',i);
%             end
            
            % copied from explor:
            
%             fn4D_dbstack
%             if isempty(N.D.selshow)
%                 delete(findobj(D.ha,'tag','ActDispIm_Sel'))
%                 return
%             end
            
            % some params
%             si = N.D.SI;
%             seldimsnum = N.D.seldims-'w';
%             selectionmarks = si.selection.getselset(seldimsnum).singleset;
%             nsel = length(selectionmarks);
            
            % display set...
            
            % before all, check that selections can be displayed
            selectionaxis = [N.D.layout.dim_locations{N.selectiondim}];
            if ~ismember(selectionaxis, {'x' 'y' 'xy' 'yx'})
                disp('selections cannot be displayed')
                delete([N.D.seldisp{ind}])
                N.D.seldisp(ind) = [];
                return
            end

            
            if nargin<2, flag = 'all'; end
            if fn_ismemberstr(flag,{'all','reset'})
                % 'findobj' allows a cleanup when some objects were not
                % removed correctly
                deleteValid(N.selectiondisplay)
                N.selectiondisplay = cell(1,length(N.selection));
                for k=1:length(N.selection)
                     displayonesel(N,k,'new');
                end
                return
            end
            
            % or display update
%             if ~isempty(N.D.curselprev) && ~isempty(strfind(N.D.selshow,'number'))
%                 set(N.D.seldisp{N.D.curselprev}(1),'color','w')
%             end
            switch flag
                case 'new'
                    for idx=ind
                        displayonesel(N,idx,'new'); 
                    end
                case {'add','change','affinity'}
                    % might be several indices
                    for k=ind, displayonesel(N,k,'pos'); end
                case 'changereferential'
                    % it is not the positions of selections that have
                    % changed, but the referential of these positions
                    % relative to the main display axes: we need then to
                    % recompute the positions of each selection display
                    % inside this new referential
                    for k = 1:length(N.selection), displayonesel(N, k, 'pos'), end
                case 'remove'
                    delete([N.D.seldisp{ind}])
                    N.D.seldisp(ind) = [];
                    nsel = length(N.D.seldisp);
                    if nsel==0, return, end
                    updateselorderdisplay(N.D)
                case 'reorder'
                    perm = value;
                    N.D.seldisp = N.D.seldisp(perm);
                    updateselorderdisplay(N.D)
                case 'indices'
                    % nothing to do
            end
%             if ~isempty(N.D.currentselection) && ~isempty(strfind(N.D.selshow,'number'))
%                 set(N.D.seldisp{N.D.currentselection}(1),'color','r')
%             end
                
        end
        function displayonesel(N,k,flag,varargin)
            % function displayonesel(D,k,'new')       - selection at index k is new
            % function displayonesel(D,k,'pos')       - has changed
            % function displayonesel(D,k,'isel',isel) - has its number changed
            % function displayonesel(D,k,'edit')      - edit mode has been turned on
            
            % flags: what kind of update
            [flagnew, flagpos, flagisel, flagedit] = ...
                fn_flags('new','pos','isel','edit',flag);
                
            %selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            %selij = selectionmarks(k);
            selij = N.selection(k);
            
            % TODO: the selection must have the information of the dimension to which
            % it applies. Maybe this has to be done with the bank to
            % synchronize filters like the zoomFilters. For now seldimsnum
            % is equal to the number of dimension information
            seldim = N.selectiondim;
            if flagnew || flagedit || flagpos
                % convert selection to displayable polygon
                selectionaxis = [N.D.layout.dim_locations{seldim}];
                if ~ismember(selectionaxis, {'x' 'y' 'xy' 'yx'})
                    error('selection cannot be displayed')
                end
                
                [polygon center] = N.D.graph.selectionMark(seldim,selij);
                %                 % visible part of the polygon
                %                 polygon = visiblePolygon(N, selij.polygon,[1, 2]);                
            end
            if flagnew || flagedit || flagisel
                colors = fn_colorset;
                col = colors(mod(k-1,size(colors,1))+1,:);
                visible = 'on';
            end
            if flagisel
                isel = varargin{1};
                str = num2str(isel);
            elseif flagnew
                str = num2str(k);
            end
            
            % Create / update objects
            if flagnew
                
                hl = [];
                %if strfind(D.selshow,'number')
                if true
                    hl(end+1) = text(center(1),center(2),str, ...
                        'Parent',N.D.ha,'color','w','visible',visible, ...
                        'horizontalalignment','center','verticalalignment','middle', ...
                        'color','w');
                    %'color',fn_switch(k==D.currentselection,'r','w'));
                end
                %if strfind(D.selshow,'shape')
                if true
                    hl(end+1) = line(polygon(1,:),polygon(2,:),'Parent',N.D.ha, ...
                        'Color',col, ...
                        'UserData',k); % set user data because this line will be used when in seledit mode
                end
                %if strfind(D.selshow,'cross')
                if false
                    hl(end+1) = line(center(1),center(2),'Parent',N.D.ha, ...
                        'Color',col,'LineStyle','none', ...
                        'Marker','+','MarkerSize',4);
                end
                set(hl,'tag','ActDispIm_Sel','HitTest','off')
                
                if k <= length(N.selectiondisplay)
                    if isgraphics(N.selectiondisplay{k})
                        delete(N.selectiondisplay{k});
                    end
                end
                
                 N.selectiondisplay{k} = hl;

                
                
            else
                hl = N.selectiondisplay{k};
                i=1; ht=[]; hs=[]; hc=[];
%                 if strfind(D.selshow,'number'), ht=hl(i); i=i+1; end
%                 if strfind(D.selshow,'shape'),  hs=hl(i); i=i+1; end
%                 if strfind(D.selshow,'cross'),  hc=hl(i); i=i+1; end
                if true, ht=hl(i); i=i+1; end
                if true,  hs=hl(i); i=i+1; end
                if false,  hc=hl(i); i=i+1; end
                he = hl(i:end);
                if flagpos
                    set(ht,'position',center)
                    set(hs,'xdata',polygon(1,:),'ydata',polygon(2,:))
                    set(hc,'xdata',center(1),'ydata',center(2))
                elseif flagisel
                    set(ht,'string',str)
                    set([hs hc he],'color',col)
                end
            end
            
            %             % Advanced selection mode (in this mode, D.seldisp = [ht hl he]
            %             % because D.selshow = 'number+shape')
            %             if ~D.seleditmode || flagisel, return, end
            %             desc = [];
            %             switch selectionmarks(k).type
            %                 case {'poly2D','mixed','point2D','line2D'} % TODO: not sure about 'point2D'
            %                     polymark = polygon;
            %                 case 'rect2D'
            %                     polymark = polygon(:,1:4); % the 5th point of polygon is a repetition of the 1st one
            %                     desc = [sel.shapes.points' sel.shapes.vectors'];
            %                 case {'ellipse2D' 'ring2D'}
            %                     c = sel.shapes.points;
            %                     u = sel.shapes.vectors;
            %                     e = sel.shapes.logic;
            %                     polymark = [c-u c+u];
            %                     desc = {c u e};
            %                 otherwise
            %                     error programming
            %             end
            %             if flagnew || flagedit
            %                 % right now, hl has 2 elements: number and shape
            %                 set(hl(2),'hittest','on','buttondownfcn', ...
            %                     @(h,evnt)seleditaction(D,get(h,'userdata'),'line'))
            %                 hl(3) = line(polymark(1,:),polymark(2,:),'Parent',D.ha, ...
            %                     'Color',col,'tag','ActDispIm_Sel', ...
            %                     'LineStyle','none','marker','.', ...
            %                     'UserData',k,'hittest','on','buttondownfcn',...
            %                     @(h,evnt)seleditaction(D,get(h,'userdata'),'point'));
            %                 if ~isempty(desc),
            %                     setappdata(hl(3),'description',desc)
            %                 end
            %                 D.seldisp{k} = hl;
            %             else
            %                 set(hl(3),'xdata',polymark(1,:),'ydata',polymark(2,:));
            %                 if ~isempty(desc)
            %                     setappdata(hl(3),'description',desc)
            %                 end
            %             end
        end
        function deletedisplayonesel(N,k)
            delete(N.selectiondisplay{k});
        end
        function output = visiblePolygon(N,points,selectionDimension)
                    % visiblePolygon(N, points);
                    
                    % @param points: nxn double list of points with
                    % coordinates vertical
                    % @param selectionDimension: 1xn list of dimension to
                    % which the selection applies
                    % @return output: same double list of points with NaN instead of
                    % points not visible and additional points to display
                    % selection properly
                    
                    % coordinates of the polygon inside a vignette (in index
                    % coordinate system)
                    %points = selij2.shapes.points;
                    np = size(points, 2);
                    
                    % replace points that are outside of vignette by NaNs
                    zoomSliceValues = N.graph.getZoom(selectionDimension,'indices');
                    
                    % set the ouput to zeros (they will be set to one if one of the
                    % dimension if it's out of display)
                    polygonIsOutOfDisplay = zeros(1,size(points,2));
                    
                    for dimension = 1:size(zoomSliceValues,2)     
                       % is equal to one if is out of limits of the zoom or if the
                       % previous value was already 1
                        polygonIsOutOfDisplay = points(dimension,:)<(zoomSliceValues(1,dimension)-.5) | points(dimension,:)>(zoomSliceValues(2,dimension)+.5) | polygonIsOutOfDisplay;

                    end
                    
                    % array of ones' boundaries in polygonIsOutOfDisplay.
                    % Represent the start (1) and finish (-1) of lines not displayed
                    boundaries=diff([0 polygonIsOutOfDisplay 0]);
                    
                    % if the first point and last point are hidden, don't
                    % consider them as start and finish of group of ones
                    if boundaries(1) == 1 && boundaries(end) == -1
                       boundaries(1)=0;
                       boundaries(end)=0;
                    end
                    
                    boundariesIndexes=find(boundaries==-1 | boundaries==1);
                    boundariesValues=nan(size(points,1),size(boundariesIndexes,2));

                    % add intermediate points between points displayed and 
                    % points not displayed
                    numberOfPointsAdded=0;
                    for boundariesIndex = 1:length(boundariesIndexes)
                        if boundariesIndexes(boundariesIndex) == 1
                            boundaryPrev = length(points)-1;
                        else
                            boundaryPrev = boundariesIndexes(boundariesIndex)-1;
                        end

                        if(boundaries(boundariesIndexes(boundariesIndex))==1)
                            % find the lowest ratio before limit of the
                            % next point
                            biggestRatio = 0;
                            for dimension = 1:size(zoomSliceValues,2)
                                vector = points(dimension,boundaryPrev) - points(dimension,boundariesIndexes(boundariesIndex));
                                vectorToLimitMin = (zoomSliceValues(1, dimension)-.5) - points(dimension,boundariesIndexes(boundariesIndex));
                                vectorToLimitMax = (zoomSliceValues(2, dimension)+.5) - points(dimension,boundariesIndexes(boundariesIndex));
                                
                                % if same sign
                                V = [vectorToLimitMax, vectorToLimitMin];
                                if (~any(diff(sign(V(V~=0)))))
                                    biggestRatio = max(biggestRatio,min(abs(vectorToLimitMin),abs(vectorToLimitMax))/abs(vector));
                                end
                            end
                            
                            boundariesValues(:,boundariesIndex) = points(:,boundariesIndexes(boundariesIndex))+(points(:,boundaryPrev)-points(:,boundariesIndexes(boundariesIndex)))*biggestRatio;
                            boundariesIndexes(boundariesIndex) = boundariesIndexes(boundariesIndex) + numberOfPointsAdded;

                        else
                            boundariesIndexes(boundariesIndex) = boundaryPrev;
                            if boundariesIndexes(boundariesIndex) == length(points)
                                boundaryNext = 2;
                            else
                                boundaryNext = boundariesIndexes(boundariesIndex)+1;
                            end

                            % find the lowest ratio before limit of the
                            % next point
                            biggestRatio = 0;
                            for dimension = 1:size(zoomSliceValues,2)
                                vector = points(dimension,boundaryNext) - points(dimension,boundariesIndexes(boundariesIndex));
                                vectorToLimitMin = (zoomSliceValues(1,dimension)-.5) - points(dimension,boundariesIndexes(boundariesIndex));
                                vectorToLimitMax = (zoomSliceValues(2,dimension)+.5) - points(dimension,boundariesIndexes(boundariesIndex));
                                
                                % if same sign
                                V = [vectorToLimitMax, vectorToLimitMin];
                                if (~any(diff(sign(V(V~=0)))))
                                    biggestRatio = max(biggestRatio,min(abs(vectorToLimitMin),abs(vectorToLimitMax))/abs(vector));
                                end
             
                            end
                            
                            boundariesValues(:,boundariesIndex) = points(:,boundariesIndexes(boundariesIndex))+(points(:,boundaryNext)-points(:,boundariesIndexes(boundariesIndex)))*biggestRatio;
                            boundariesIndexes(boundariesIndex) = boundaryNext + numberOfPointsAdded;
%                             points(selectionDimension,boundariesIndexes(boundariesIndex)) = points(:,boundariesIndexes(boundariesIndex))+(points(:,boundaryNext)-points(:,boundariesIndexes(boundariesIndex)))*biggestRatio;
                        end
                        
                        numberOfPointsAdded = numberOfPointsAdded + 1;
                    end

                    % if column n of polygonIsOutOfDisplay is equal to 1 it
                    % means the point in column n of polygon is out of
                    % display so replace all values on those columns by Nan
                    % to don't be drawn

                    % coordinates in the full nd index coordinate system

                    points(:,polygonIsOutOfDisplay) = NaN;

                    output = ones(N.D.nd, np + size(boundariesValues,2));
                    
                    output(selectionDimension,setdiff(1:end,boundariesIndexes)) = points;
                    output(selectionDimension,boundariesIndexes) = boundariesValues;

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
            N.addListener(Z,'ChangedZoom',@(u,e)set(obj,'value',Z.zoomvalue));
        end
    end
    
    methods
       
        %  if one of the dimension of the cross is hidden, hide the
        % cross center as well
        function updateCrossCenterVisibility(N)
            if fn_switch(get(N.cross(1), 'Visible')) && fn_switch(get(N.cross(2), 'Visible'))
                set(N.cross(3), 'Visible', 'on')
            else
                set(N.cross(3), 'Visible', 'off')
            end
        end
        
        
        % return true if point (graph coordinates) is part of the slice
        % data displayed by converting the point to slice coordinates and
        % test if its between minimal and maximal values for all dimensions
        %
        % @param point: 2xn double 
        % @param perdim: whether to return a single logical value per point
        %                (point is/isn't inside slice, default behavior)
        %                or a vector (for each dimension, if perdim=true)
        % @return output: 1xn boolean
        function output = isOutOfDisplay(N, point, perdim)
            if nargin<3, perdim = false; end
            
            % get slice values of the point
            ijk = N.graph.graph2slice(point,'invertible',true)'; % n x ndim
            
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