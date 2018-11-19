classdef displaygraph < handle
    
    properties (SetAccess='private')
        % parent xplr.viewdisplay object
        D
        % graphic objects
        xyticks
    end
    properties (Dependent, SetAccess='private')
        ha
        org
    end
    % parameters
    properties (SetAccess='private')
        xsep = .2;
        ysep = .2;
    end
    % pre-computations
    properties (SetAccess='private')
        zslicesz    % current size of the zoomed slice
        filling     % how much of the space available for each dimension if filled (vector of values <=1)
    end
    properties (Access='private')
        steps
    end
        
    
    % Constructor, Get dependent
    methods
        function G = displaygraph(D)
            % parent xplr.viewdisplay object
            G.D = D;
        end
        function ha = get.ha(G)
            ha = G.D.ha;
        end
        function org = get.org(G)
            org = G.D.org;
        end
    end
    
    % Pre-computations
    methods (Access='private')
        function [st sz fill xpair ypair] = computeStepsPrivate(G,orgin)
            % actual steps computation occurs here
            % size and organization
            dosignal = strcmp(G.D.displaymode,'time courses'); % slight difference for x between time courses and images
            header = G.D.zslice.header;
            sz = [header.n]; nd = length(sz);
            szorig = G.D.slice.sz; okdim = (szorig>1);
            xorg = orgin.x;         yorg = orgin.y;
            nx = length(xorg);      ny = length(yorg);
            zoom = G.getZoom('value');              % the 'true' coordinates (i.e. in data before zooming) of space edges
            [idxoffset bin] = G.getZoom('off&bin'); % coordinates conversions between original and zoomed data
            
            % zooming
            % It shall be noted that while the zoom specification are
            % arbitrary real values, the data will be cut at integer
            % values, therefore there are two ways to display it:
            % - either shifting it by the mismatch between these real /
            %   integer values (in this case moving the zoom will result in
            %   smooth movements)
            % - either occupying the whole dedicated space without caring
            %   about this mismatch (in this case moving the zoom will
            %   result in step movements)
            % The most pleasant navigation in the data seems to be obtained
            % when the first choice is made for the "most outside"
            % dimensions, and the second is made for the "inner"
            % dimensions.
            % (convert real zoom values into zoomed data coordinates)
            % we have: xtrue   = .5 + offset + (xzoomed-.5)*bin
            % and:     xzoomed = .5 + (xtrue-.5-offset)/bin
            zm = .5 + (mean(zoom)-.5-idxoffset)./bin;   % zoom centers in zoomed data coordinates
            ze = diff(zoom)./bin;                       % zoom extents in zoomed data coordinates
            % (these real zoom values will not be used for the "internal"
            % coordinates; they will neither be used for the "2D grid"
            % organization, but this does not appear in the lines below)
            if any(orgin.xy) || any(orgin.yx)
                din = [xorg yorg];
            else
                lastx = find(okdim(xorg),1,'last'); % index of "extern" x dimension
                lasty = find(okdim(yorg),1,'last'); % index of "extern" y dimension
                din = [xorg(1:lastx-1) yorg(1:lasty-1)]; % "internal" dimensions
            end
            zm(din) = (sz(din)+1)/2;
            ze(din) = sz(din);
            % (in addition, a small correction is applied for signals,
            % whose edges will be at 1 and n instead of .5 and n+.5 for
            % images and grid cells)
            if dosignal && nx>=1 && sz(xorg(1))>1
                ze(xorg(1)) = ze(xorg(1))-1;
            end 

            % specific data aspect ratio for some dimensions?
            xhead = header(xorg);   yhead = header(yorg);
            [xpair ypair] = checkpairs(xhead,yhead,dosignal);
            [w h] = fn_pixelsize(G.D.ha);
            axisratio = h/w;

            % G.filling will get values assigned both for 'grid' and
            % 'linear' arrangements
            [fill st.xspan st.yspan] = deal(zeros(1,nd),zeros(1,nx),zeros(1,ny));
            
            % does one dimension have 2D grid organization
            xavail = 1; yavail = 1; % total available x- and y-span
            if any(orgin.xy) || any(orgin.yx)
                % which dimension
                if any(orgin.xy)
                    st.xydim = orgin.xy;
                    xymode = 'xy';
                else
                    st.xydim = orgin.yx;
                    xymode = 'yx';
                end
                
                % TODO: the lines below are not correct, remove?
                % % span
                % [st.xspan(st.xydim) st.yspan(st.xydim)] = deal(xavail,yavail);
                                
                % determine number of column: what aspect ratio is desired
                % for the grid elements?
                nelem = header(st.xydim).n;
                if nx && ny && xpair(end) && ypair(end)
                    elemratio = (yhead(end).scale*yhead(end).n)/(xhead(end).scale*xhead(end).n);
                else
                    elemratio = 1; % this value will be only loosely respected
                end
                if strcmp(xymode,'xy')
                    ncol = fn_coerce(round(sqrt(nelem*elemratio*xavail/yavail/axisratio)),[1 nelem]);
                    nrow = ceil(nelem/ncol);
                else
                    nrow = fn_coerce(round(sqrt(nelem/elemratio/xavail*yavail*axisratio)),[1 nelem]);
                    ncol = ceil(nelem/nrow);
                end
                
                % set grid positions
                st.xyoffsets = zeros(2,nelem);
                i0 = (ncol+1)/2; x_step = xavail/ncol;
                j0 = (nrow+1)/2; y_step = -yavail/nrow;
                st.xysteps = [x_step y_step];
                for k=1:nelem
                    switch xymode
                        case 'xy'
                            [i j] = ind2sub([ncol nrow],k);
                        case 'yx'
                            [j i] = ind2sub([nrow ncol],k);
                    end
                    st.xyoffsets(:,k) = [(i-i0)*x_step; (j-j0)*y_step];
                end
                if szorig(st.xydim)>1
                    xavail = xavail/ncol/(1+G.xsep);
                    yavail = yavail/nrow/(1+G.ysep);
                end
                fill(st.xydim) = 1;
            else
                st.xydim = [];
                st.xyoffsets = zeros(2,1);
            end
            
            % define steps while ensuring aspect ratio for pairs (start
            % from the last dimensions)
            [st.xoffset st.xstep] = deal(zeros(1,nx)); % offsets will be adjusted so as to keep data point of "center-zoom" coordinates in the middle of the display
            [st.yoffset st.ystep] = deal(zeros(1,ny));
            fill([xorg yorg]) = 1;
            ix = nx; iy = ny;
            while ix>0 || iy>0
                % go down to the next pair
                ixnext = find(xpair(1:ix),1,'last');
                iynext = find(ypair(1:iy),1,'last');
                if isempty(ixnext), [ixnext iynext] = deal(1); end
                % (x)
                for ix=ix:-1:ixnext
                    d = xorg(ix);
                    st.xspan(ix) = xavail;
                    st.xstep(ix) = xavail / ze(d);
                    st.xoffset(ix) = -zm(d)*st.xstep(ix);   % middle of zoom should be placed at the middle of the available space
                    if szorig(d)>1, xavail = st.xstep(ix) / (1+G.xsep); end  % available x-span for (ix-1)th dimension
                end
                % (y)
                for iy=iy:-1:iynext
                    d = yorg(iy);
                    st.yspan(iy) = yavail;
                    st.ystep(iy) = -yavail / ze(d);        % start from top of the screen (i.e. higher values of y) rather than bottom
                    st.yoffset(iy) = -zm(d)*st.ystep(iy);   % middle of zoom should be placed at the middle of the available space
                    if szorig(d)>1, yavail = abs(st.ystep(iy)) / (1+G.ysep); end % available y-span for (iy-1)th dimension
                end
                
                % arrange values to maintain aspect ratio for the pair if
                % there is a pair
                if isempty(ix) || (ix==1 && ~xpair(ix)), break, end
                curratio = abs(st.ystep(iy))/st.xstep(ix) * axisratio;
                targetratio = yhead(iy).scale/xhead(ix).scale;
                correction = targetratio/curratio;
                if correction>1
                    % need to reduce x-span
                    d = xorg(ix);
                    st.xspan(ix) = st.xspan(ix)/correction;
                    st.xoffset(ix) = st.xoffset(ix) + zm(d)*st.xstep(ix)*(1-1/correction);
                    st.xstep(ix) = st.xstep(ix)/correction;
                    fill(d) = 1/correction;
                    xavail = xavail/correction;
                elseif correction<1
                    % need to reduce y-span
                    d = yorg(iy);
                    st.yspan(iy) = st.yspan(iy)*correction;
                    st.yoffset(iy) = st.yoffset(iy) + zm(d)*st.ystep(iy)*(1-1*correction);
                    st.ystep(iy) = st.ystep(iy)*correction;
                    fill(d) = correction;
                    yavail = yavail*correction;
                end
                ix = ix-1;
                iy = iy-1;
            end
            
            % case empty
            if isempty(xorg)
                % coordinate 1 should go to the center
                st.xoffset = -xavail;
                st.xstep = xavail;
            end
            if isempty(yorg)
                % coordinate 1 should go to the center
                st.yoffset = yavail;
                st.ystep = -yavail; % ystep must be negative
            end
        end
    end
    methods
        function anychg = computeSteps(G)
            % function anychg = computeSteps(G)
            %---
            % sets G properties xoffset xstep yoffset ystep and tells
            % whether they were changed
            
            % compute steps
            if nargout>0, prevsteps = G.steps; end
            [G.steps G.zslicesz G.filling xpair] = computeStepsPrivate(G,G.org); %#ok<ASGLU>
            
            % any change
            if nargout>0
                 anychg = ~isequal(G.steps,prevsteps);
            end
        end
    end
    
    % Ticks
    methods
        function setTicks(G)
            axsiz = fn_pixelsize(G.ha);
            axsizinch = axsiz/get(0,'ScreenPixelsPerInch');
            targetspacinginch = .5; % optimal space between ticks in inches
            maxnarrow = 2; % maximal further narrowing of this optimal space
            st = G.steps;
            
            % remove previous xy ticks
            deleteValid(G.xyticks)
            G.xyticks = [];
            
            % stop if data is too large for being displayed
            if G.D.nodisplay, return, end
            
            % x and y
            for k = 1:2
                switch k
                    case 1
                        [f ff] = deal('x','yx');
                    case 2
                        [f ff] = deal('y','xy');
                end
                d = G.D.activedim.(f);
                if isempty(d)
                    tick = [];
                    ticklabels = {};
                else
                    head = G.D.zslice.header(d);
                    n = head.n;
                    % conversion between data coordinates and graph
                    domeasure = head.ismeasure;
                    if d==G.org.(ff)
                        % step for ncol (or nrow) data points
                        ncol = find(diff(st.xyoffsets(k,:)),1);
                        if isempty(ncol), ncol = size(st.xyoffsets,2); end
                        f_step = st.xysteps(k); 
                        % step for one data point
                        f_step = f_step / ncol;
                        f_off = st.xyoffsets(k,1) - (ncol+1)/2*f_step;
                        domeasure = false;
                    else
                        jf = find(d==G.org.(f),1);
                        if isempty(jf), error('%s activedim must have either ''%s'' or ''%s'' organization!',f,f,ff), end
                        switch f
                            case 'x'
                                f_off = st.xyoffsets(k,1) + st.xoffset(jf) + sum(st.xoffset(jf+1:end)+st.xstep(jf+1:end));
                                f_step = st.xstep(jf);
                            case 'y'
                                f_off = st.xyoffsets(k,1) + st.yoffset(jf) + sum(st.yoffset(jf+1:end)+st.ystep(jf+1:end));
                                f_step = st.ystep(jf);
                        end
                    end
                    % measure or labels?
                    if domeasure
                        % target space between ticks
                        targetspacing = targetspacinginch / axsizinch(k);   % target spacing in axes coordinates
                        fspan = fn_switch(f,'x',st.xspan,'y',st.yspan);
                        targetspacing = targetspacing/min(1/fspan(jf),maxnarrow); % let this target increase up to a factor of two when dimension occupies only a fraction of the space
                        % target space in data coordinates
                        [start scale] = deal(head.start,head.scale);
                        target = targetspacing / abs(f_step) * scale;
                        % actual step that will be used
                        t10 = log10(target);
                        tests = [1 2 5 10];
                        [~, idx] = min(abs(mod(t10,1)-log10(tests)));
                        step = 10^floor(t10) * tests(idx);
                        % tick values
                        ticksdata = step * (ceil((start-.5*scale)/step):floor((start+(n-.5)*scale)/step)); % data coordinates
                        if f=='y', ticksdata = fliplr(ticksdata); end % make it ascending order
                        ticksidx = 1 + (ticksdata-start)/scale; % data indices coordinates
                        % tick labels
                        ticklabels = fn_num2str(ticksdata,'cell');
                    else
                        % ticks for each data point
                        ticksidx = 1:n;
                        ticklabels = row(head.getItemNames());
                        % 2D grid: put text rather than using axes ticks
                        % system
                        if d==G.org.(ff)
                            xy = fn_add([0; -st.xysteps(2)/2],st.xyoffsets);
                            G.xyticks = gobjects(1,n);
                            for i=1:n
                                G.xyticks(i) = text(xy(1,i),xy(2,i),ticklabels{i}, ...
                                    'parent',G.ha,'hittest','off', ...
                                    'horizontalalignment','center'); %,'verticalalignment','top')
                            end
                            ticksidx = []; ticklabels = {};
                        end
                        if f=='y'
                            ticksidx = fliplr(ticksidx); 
                            ticklabels = fliplr(ticklabels);
                        end
                    end
                    tick = f_off + ticksidx*f_step;
                end
                % set ticks!
                if strcmp(f,'x')
                    set(G.ha,'xtick',tick,'xticklabel',ticklabels)
                else
                    set(G.ha,'ytick',tick,'yticklabel',ticklabels)
                end
            end
            
        end
    end

    % Coordinates conversions
    methods
        function [zoom bin] = getZoom(G,varargin)
            % function zoom = getZoom(G[,dim][,'value|effective|indices'])
            % function [offset bin] = getZoom(G[,dim],'off&bin')
            %---
            % Get zoom value for specified dimensions. 
            % Different modes affect the returned value:
            % - 'value' [default]   returns the zooming value specified to
            %               the zoom filter
            % - 'effective' zooming value, after taking into account the
            %               extra space due to not completely filling the
            %               available space (to preserve data aspect ratio)
            % - 'indices'   returns zoom values of the zoomed data points
            % - 'off&bin'   returns offset (i.e. number of data points from
            %               the beginning which are not shown) and binning
            %               value
            
            % input
            dim = 1:G.D.nd;
            mode = 'value';
            for i=1:length(varargin)
                a = varargin{i};
                if isnumeric(a)
                    dim = a;
                else
                    mode = a;
                end
            end
            
            % output
            zfilters = G.D.zoomfilters(dim);
            switch mode
                case {'value' 'effective'}
                    zoom = cat(1,zfilters.zoomvalue)';
                    if strcmp(mode,'effective')
                        zoom = fn_add(mean(zoom), fn_mult([-.5; .5],diff(zoom)./G.filling(dim)));
                    end
                case 'indices'
                    zoom = zeros(2,length(dim));
                    for i=1:length(dim), zoom(:,i) = zfilters(i).indicesout([1 end]); end
                case 'off&bin'
                    offset = zeros(1,length(dim));
                    for i=1:length(dim), offset(i) = zfilters(i).indicesin(1)-1; end
                    zoom = offset;
                    bin = [zfilters.bin];
            end
        end
        function xy = zslice2graph(G,ijk)
            % function xy = zslice2graph(G,ijk)
            %---
            % convert data coordinates to graph positions (btw 0 and 1)
            % Input ijk can be a vector (coordinates of a single point) or
            % array (coordinates of as many points as columns).
            st = G.steps;
            
            if isvector(ijk), ijk = ijk(:); end
            x = sum(fn_add(st.xoffset(:), fn_mult(ijk(G.org.x,:),st.xstep(:))),1);
            y = sum(fn_add(st.yoffset(:), fn_mult(ijk(G.org.y,:),st.ystep(:))),1);
            xy = [x; y];
            if ~isempty(st.xydim), xy = xy + st.xyoffsets(:,ijk(st.xydim,:)); end
        end
        function ijk = graph2zslice(G,xy)
            % function ijk = graph2zslice(G[,xy])
            if nargin<2
                xy = get(G.D.ha,'currentpoint');
                xy = xy(1,1:2);
            end
            if isvector(xy), xy = xy(:); end
            np = size(xy,2);
            ijk = zeros(G.D.nd,np);
            st = G.steps;
            sz = G.zslicesz;
            
            % xy/yx
            if ~isempty(st.xydim)
                % take advantage on the fact that the grid spans the full
                % axis
                d = st.xydim;
                x = .5+xy(1,:); y = .5-xy(2,:);
                ncol = round(1/st.xysteps(1)); icol = .5 + x*ncol;
                nrow = round(1/abs(st.xysteps(2))); irow = .5 + y*nrow;
                if d==G.org.xy
                    ijk(d,:) = icol + ncol*round(irow-1);
                else
                    ijk(d,:) = irow + nrow*round(icol-1);
                end
                xy = xy - st.xyoffsets(:,fn_coerce(round(ijk(d,:)),1,sz(d)));
            end
            
            % x 
            x = xy(1,:);
            xorg = G.org.x;
            for ix = length(xorg):-1:1
                d = xorg(ix);
                x = x - st.xoffset(ix);
                ijk(d,:) = fn_div(x,st.xstep(ix));
                x = x - fn_mult(round(ijk(d,:)),st.xstep(ix));
            end
            
            % y
            y = xy(2,:);
            yorg = G.org.y;
            for iy = length(yorg):-1:1
                d = yorg(iy);
                y = y - st.yoffset(iy);
                ijk(d,:) = fn_div(y,st.ystep(iy));
                y = y - fn_mult(round(ijk(d,:)),st.ystep(iy));
            end
        end
        function ijk = graph2slice(G,xy)
            % function ijk = graph2slice(G[,xy])
            if nargin<2
                xy = get(G.D.ha,'currentpoint');
                xy = xy(1,1:2);
            end
            
            % coordinates in zoomed slice
            zijk = graph2zslice(G,xy);
            
            % convert to before zooming
            [idxoffset bin] = G.getZoom('off&bin');
            ijk = fn_add(idxoffset(:), .5+fn_mult(zijk-.5,bin(:)));
        end
        function M = gettransform(G,ijk,ylim_or_ybase,yextent)
            % function M = gettransform(G,ijk,ybase)
            %---
            % matrix transformation to place curve/image at ijk data
            % coordinates; note that only coordinates not belonging to the
            % curve/image will be taken into account.
            % 
            % Input:
            % - ijk     must be a vector (coordinates of a single point)
            % - ylim    (for 'time courses' mode only) 2 data values
            %           corresponding respectively to the bottom and top
            %           positions in the space dedicated to displaying this
            %           curve
            % - ybase, yextent  alternate way of specifying ylim (ylim =
            %           ybase + [-.5 .5]*ystep)
            st = G.steps;

            % Initialize matrix
            M = repmat(eye(4),[1 1 size(ijk,2)]);
            
            % Scale: depends only on in-curve/in-image dimension(s)
            % (x)
            xscale = st.xstep(1);
            M(1,1,:) = xscale;
            % (y)
            switch G.D.displaymode
                case 'image'
                    % not possible to have negative values -> orienting the
                    % images downward will be achieved by inverting y
                    % coordinates at the stage of patch creation
                    yscale = abs(st.ystep(1));
                case 'time courses'
                    if nargin==3
                        ylim = ylim_or_ybase;
                        ybase = mean(ylim);
                        yextent = diff(ylim);
                    else
                        ybase = ylim_or_ybase;
                    end
                    yscale = abs(st.ystep(1) / yextent); % y-values of time courses should not be oriented downward
                otherwise
                    error 'invalid display mode'
            end
            M(2,2,:) = yscale; 
            
            % Offset: handle separately offsets relative to
            % in-curve/in-image  dimension(s) and to other dimensions
            % (x)
            M(1,4,:) = fn_add( sum(st.xoffset), sum(fn_mult(column(st.xstep(2:end)),ijk(G.org.x(2:end),:)),1) );
            % (y)
            switch G.D.displaymode
                case 'image'
                    M(2,4,:) = fn_add( sum(st.yoffset), sum(fn_mult(column(st.ystep(2:end)),ijk(G.org.y(2:end),:)),1) );
                case 'time courses'
                    % we subtract yscale*ybase so that ybase goes to
                    % the center
                    M(2,4,:) = fn_add( sum(st.yoffset)-yscale*ybase, sum(fn_mult(st.ystep(:),ijk(G.org.y,:)),1) );
            end
            % (xy)
            if ~isempty(st.xydim)
                M(1:2,4,:) = M(1:2,4,:) + permute(st.xyoffsets(:,ijk(st.xydim,:)),[1 3 2]); 
            end
        end
        function pos = labelPosition(G,dim,orgin)
            % function pos = labelPosition(G,d[,orgin])
        
            % steps
            if nargin==3
                st = computeStepsPrivate(G,orgin);
            else
                orgin = G.org;
                st = G.steps;
                if isempty(st)
                    % code here added on 04/06/2018 to avoid error
                    if G.D.nodisplay
                        warning 'Please edit code to better handle this case'
                        pos = zeros(1,length(dim));
                        return
                    else
                        error 'programming'
                    end
                end
            end
            
            % label positions
            n = length(dim);
            pos = zeros(1,n);
            for i=1:n
                d = dim(i);
                if ~isempty(orgin.x) && d==orgin.x(end)
                    pos(i) = st.xyoffsets(1,1);
                elseif any(d==orgin.x)
                    ix = find(d==orgin.x,1);
                    pos(i) = st.xyoffsets(1,1) + sum(st.xoffset(ix+1:end) + st.xstep(ix+1:end));
                    if pos(i)<-.5, pos(i) = pos(i) + sum(st.xstep(ix+1:end)); end % first grid element is more than half-outside
                elseif ~isempty(orgin.y) && d==orgin.y(end)
                    pos(i) = st.xyoffsets(2,1);
                elseif any(d==orgin.y)
                    iy = find(d==orgin.y,1);
                    pos(i) = st.xyoffsets(2,1) + sum(st.yoffset(iy+1:end) + st.ystep(iy+1:end));
                    if pos(i)>.5, pos(i) = pos(i) + sum(st.ystep(iy+1:end)); end % first grid element is more than half-outside
                elseif d==st.xydim
                    pos(i) = 0;
                else
                    pos(i) = 0;
                end
            end
        end
    end
end



%---
function [xpair ypair] = checkpairs(xhead,yhead,dosignal)
% this function is very ad-hoc and could be improved

nx = length(xhead);  ny = length(yhead);
xpair = zeros(1,nx); ypair = zeros(1,ny);

% restrict to dimensions which are measures
xok = false(1,nx);   yok = false(1,ny);
for i=1:nx, xok(i) = xhead(i).ismeasure; end
for i=1:ny, yok(i) = yhead(i).ismeasure; end
if ~any(xok) || ~any(yok), return, end

% look for dimensions being in the same space!
xunits = {xhead.unit};
yunits = {yhead.unit};
if xok(1) && yok(1) && isequal(xunits(1),yunits(1))
    [xpair(1) ypair(1)] = deal(1,1);
    [xok(1) yok(1)] = deal(false);
end
if nx>=2 && ny>=2 && xok(nx) && yok(ny) && isequal(xunits(nx),yunits(ny))
    [xpair(nx) ypair(ny)] = deal(ny,nx);
    [xok(nx) yok(ny)] = deal(false);
end
if dosignal && nx>=2 && xok(2) && yok(1) && isequal(xunits(2),yunits(1))
    % (x1 cannot be paired with y2 for images)
    [xpair(2) ypair(1)] = deal(1,2);
    [xok(1) yok(1)] = deal(false);
end
if nx>=2 && xok(2) && ny>=2 && yok(2) && isequal(xunits(2),yunits(2))
    [xpair(2) ypair(2)] = deal(2,2);
    [xok(2) yok(2)] = deal(false); %#ok<NASGU>
end

end
