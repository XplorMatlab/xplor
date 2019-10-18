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
            dim = [N.D.activedim.x N.D.activedim.y];
            if isempty(dim), return, end
            switch get(N.hf,'SelectionType')
                case 'normal'
                    % zoom in or select point
                    if pointonly                        
                        dozoom = false;
                    else
                        rect = fn_mouse(N.ha,'rectangle-');
                        dozoom = any(any(diff(rect,1,2)));
                    end
                    if dozoom
                        ijk = N.graph.graph2slice(rect(:,[1 3]));
                        zoom = ijk(dim,:)';
                        for i=1:length(dim), zoom(:,i) = sort(zoom(:,i)); end
                        N.D.zoomslicer.setZoom(dim,zoom)
                    else
                        N.crossCenter = point;
                    end
                case 'open'
                    % zoom reset
                    zoom = repmat(':',1,length(dim));
                    N.D.zoomslicer.setZoom(dim,zoom)
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
                N.addListener(P,'ChangedPoint',@(u,e)movedPoint(N, d))
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
        function movedPoint(N, d)
            F = N.dimfilters{d};
            if isempty(F)
                error('movedPoint callback called for dimension %i, but there is no filter in this dimension!', d)
            end
            P = F.P;
            ijk = N.graph.graph2slice(N.crossCenter);
            ijk(d) = P.index;
            N.crossCenter = N.graph.slice2graph(ijk);
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
            
            % update the point filters
            try
                ijk = N.graph.graph2slice(crossCenter,true);
                for d = 1:length(ijk)
                    P = N.dimfilters{d}.P;
                    if ~isempty(P)
                        P.index = ijk(d);
                    end
                end
            end
            
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
                        N.crossCenter(1) = p(1);
                    case 2
                        N.crossCenter(2) = p(2);
                    case 3
                        N.crossCenter = p;
                    otherwise
                        error('wrong il')
                end
                
                %if do1d
                %if il~=1
                
                %end
                %if il~=2
                %si.ij2 = AX2IJ(si,p(1));
                %end
                %else
                %   ij2 = AX2IJ(si,p([1 2])');
                %   switch il
                %       case 1 % move x only
                %           si.ij2(1) = ij2(1);
                %       case 2 % move y only
                %            si.ij2(2) = ij2(2);
                %        case 3 % move x and y
                %           si.ij2 = ij2;
                %   end
                %end
            end
            
        end
        function removecross(N)
            delete(N.cross)
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
            origin = row(N.graph.graph2slice()); % current point in data coordinates
            zoomfactor = 1.5^nscroll;
            dim = [N.D.activedim.x N.D.activedim.y];
            if nscroll<0 && ~any([N.D.org.xy N.D.org.yx])
                % it does not make sense to zoom-in in a dimensions which
                % does not fill its available space due to aspect ratio
                % constraints
                dim(N.graph.filling(dim)<1) = [];
            end
            if isempty(dim), return, end
            zoom = N.graph.getZoom(dim); %,'effective');
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
    
end