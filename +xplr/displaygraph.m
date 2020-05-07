classdef displaygraph < xplr.graphnode
    % This class takes care of conversions between slice and graph coordinates,
	% and of x/y ticks

    % properties
    properties (SetAccess='private')
        D % parent xplr.viewdisplay object
        % graphic objects
        xyticks 
        separationlines
    end
    properties (SetObservable, AbortSet=true)
        showseparation = false;
        separationcolor = [.8 .8 .8];
    end
    % shortcuts
    properties (Dependent, SetAccess='private')
        ha
    end
    % parameters
    properties (SetAccess='private')
        xsep = .2;
        ysep = .2;
    end
    % pre-computations
    properties (SetAccess='private')
        layout
        zslicesz    % current size of the zoomed slice
        filling     % how much of the space available for each dimension if filled (vector of values <=1)
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
            xlayout = orgin.x;         ylayout = orgin.y;
            nx = length(xlayout);      ny = length(ylayout);
            zoom = G.getZoom('value');              % the 'true' coordinates (i.e. in data before zooming) of space edges
            [idxoffset bin] = G.getZoom('off&bin'); % coordinates conversions between original and zoomed data
            
            % zooming
            % It shall be noted that while the zoom specification are
            % arbitrary real values, the data will be cut at integer
            % values, therefore there are two ways to display it:
            % 1 - either shifting it by the mismatch between these real /
            %   integer values (in this case moving the zoom will result in
            %   smooth movements)
            % 2 - either occupying the whole dedicated space without caring
            %   about this mismatch (in this case moving the zoom will
            %   result in step movements)
           
            % Zoom method 1 "continuous":
            % (convert real zoom values into zoomed data coordinates)
            % we have: xtrue   = .5 + offset + (xzoomed-.5)*bin
            % and:     xzoomed = .5 + (xtrue-.5-offset)/bin
            zm1 = .5 + (mean(zoom)-.5-idxoffset)./bin;   % zoom centers in zoomed data coordinates
            ze1 = diff(zoom)./bin;                       % zoom extents in zoomed data coordinates
            
            % Zoom method 2 "steps":
            zm2 = (sz+1)/2;
            ze2 = sz;
            
            %             % Choose method 2 for internal coordinates anf for "2D grid"
            %             % organization. Choose method 1 otherwise.
            %             [zm ze] = deal(zm1, ze1);
            %             if any(orgin.xy) || any(orgin.yx)
            %                 din = [xlayout ylayout];
            %             else
            %                 lastx = find(okdim(xlayout),1,'last'); % index of "extern" x dimension
            %                 lasty = find(okdim(ylayout),1,'last'); % index of "extern" y dimension
            %                 din = [xlayout(1:lastx-1) ylayout(1:lasty-1)]; % "internal" dimensions
            %             end
            %             zm(din) = zm2(din);
            %             ze(din) = ze2(din);
            
            % Always choose method 2
            [zm ze] = deal(zm2, ze2);
            
            % (in addition, a small correction is applied for signals,
            % whose edges will be at 1 and n instead of .5 and n+.5 for
            % images and grid cells)
            if dosignal && nx>=1 && sz(xlayout(1))>1
                ze(xlayout(1)) = ze(xlayout(1))-1;
            end 

            % specific data aspect ratio for some dimensions?
            xhead = header(xlayout);   yhead = header(ylayout);
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
                    elemratio = abs((yhead(end).scale*yhead(end).n)/(xhead(end).scale*xhead(end).n));
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
                [st.xyncol, st.xynrow] = deal(ncol,nrow);
                
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
            fill([xlayout ylayout]) = 1;
            ix = nx; iy = ny;
            while ix>0 || iy>0
                % go down to the next pair
                ixnext = find(xpair(1:ix),1,'last');
                iynext = find(ypair(1:iy),1,'last');
                if isempty(ixnext), [ixnext iynext] = deal(1); end
                % (x)
                for ix=ix:-1:ixnext
                    d = xlayout(ix);
                    st.xspan(ix) = xavail;
                    st.xstep(ix) = xavail / ze(d);
                    st.xoffset(ix) = -zm(d)*st.xstep(ix);   % middle of zoom should be placed at the middle of the available space
                    if szorig(d)>1, xavail = st.xstep(ix) / (1+G.xsep); end  % available x-span for (ix-1)th dimension
                end
                % (y)
                for iy=iy:-1:iynext
                    d = ylayout(iy);
                    st.yspan(iy) = yavail;
                    st.ystep(iy) = -yavail / ze(d);        % start from top of the screen (i.e. higher values of y) rather than bottom
                    st.yoffset(iy) = -zm(d)*st.ystep(iy);   % middle of zoom should be placed at the middle of the available space
                    if szorig(d)>1, yavail = abs(st.ystep(iy)) / (1+G.ysep); end % available y-span for (iy-1)th dimension
                end
                
                % arrange values to maintain aspect ratio for the pair if
                % there is a pair
                if isempty(ix) || (ix==1 && ~xpair(ix)), break, end
                curratio = abs(st.ystep(iy))/st.xstep(ix) * axisratio;
                targetratio = abs(yhead(iy).scale/xhead(ix).scale);
                correction = targetratio/curratio;
                if correction>1
                    % need to reduce x-span
                    d = xlayout(ix);
                    st.xspan(ix) = st.xspan(ix)/correction;
                    st.xoffset(ix) = st.xoffset(ix) + zm(d)*st.xstep(ix)*(1-1/correction);
                    st.xstep(ix) = st.xstep(ix)/correction;
                    fill(d) = 1/correction;
                    xavail = xavail/correction;
                elseif correction<1
                    % need to reduce y-span
                    d = ylayout(iy);
                    st.yspan(iy) = st.yspan(iy)*correction;
                    st.yoffset(iy) = st.yoffset(iy) + zm(d)*st.ystep(iy)*(1-1*correction);
                    st.ystep(iy) = st.ystep(iy)*correction;
                    fill(d) = correction;
                    yavail = yavail*correction;
                end
                ix = ix-1;
                iy = iy-1;
            end
            
            % available space inside the most interior dimensions (this
            % will be used in two cases: 1) when there are no data
            % dimensions in the x- or y-axis, i.e. available space is 1;
            % 2) for time courses display, as "data value" becomes as a new
            % dimension below the most interior dimension in the y-axis.
            st.xavailable = xavail;
            st.yavailable = yavail;
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
            G.layout = G.D.layout;
            [G.steps, G.zslicesz, G.filling, xpair] = computeStepsPrivate(G,G.layout); 
            
            % any change
            if nargout>0
                 anychg = ~isequal(G.steps,prevsteps);
            end
        end
    end
    
    % Ticks
    methods (Access='private')
        function minimumspacing = ticksMinimumSpacing(G,k)
            % k = 1 for 'x', 2 for 'y'
            axsiz = fn_pixelsize(G.ha);
            axsizinch = axsiz/get(0,'ScreenPixelsPerInch');
            minimumspacinginch = .2; % optimal space between ticks in inches
            
            % target space between ticks
            minimumspacing = minimumspacinginch / axsizinch(k);   % target spacing in axes coordinates
        end
        function step = nicestep(G,minstep)
            t10 = log10(abs(minstep));
            tests = [1 2 5 10];
            [~, idx] = find(log10(tests)>=mod(t10,1),1,'first');
            step = sign(minstep) * 10^floor(t10) * tests(idx);
        end
        function [tickvalues, ticklabels] = nicevalues(G,valuestart,valuestop,minsubstep)
            % there will be two sub-steps per step (i.e. one value label
            % every two ticks)
            nsubstep = 2;
            minstep = minsubstep * nsubstep;
            % actual step that will be used: minimal "nice step" that is
            % larger than minstep
            step = G.nicestep(minstep);
            % tick values for all substeps
            substep = step/nsubstep;
            tickvalues = substep * (ceil(valuestart/substep):floor(valuestop/substep)); % data coordinates
            % tick labels only for steps
            ntick = length(tickvalues);
            ticklabels = cell(1,ntick);
            dolabel = (mod(tickvalues,step)==0);
            ticklabels(dolabel) = fn_num2str(tickvalues(dolabel),'cell');
        end
    end
    methods
        function setTicks(G)
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
                        [f, ff] = deal('x','yx');
                    case 2
                        [f, ff] = deal('y','xy');
                end
                d = G.D.activedim.(f);
                
                % no active dimension -> no tick
                if isempty(d) || G.D.V.slice.header(d).n == 1
                    switch f
                        case 'x'
                            set(G.D.ha,'xtick',[])
                        case 'y'
                            % note that method setValueTicks, which is
                            % normally called AFTER setTicks, might set
                            % yticks later
                            set(G.D.ha,'ytick',[])
                    end
                    continue
                end
                
                % header
                head = G.D.zslice.header(d);
                n = head.n;

                % conversion between data coordinates and graph
                domeasure = head.ismeasure;
                dogrid = any(d==G.layout.(ff));
                if dogrid
                    % step for ncol (or nrow) data points
                    ncol = find(diff(st.xyoffsets(k,:)),1);
                    if isempty(ncol), ncol = size(st.xyoffsets,2); end
                    f_step = st.xysteps(k); 
                    % step for one data point
                    f_step = f_step / ncol;
                    f_off = st.xyoffsets(k,1) - (ncol+1)/2*f_step;
                    domeasure = false;
                    
                    % no constraint on a minimumstep
                    minimumstep = 0;
                else
                    jf = find(d==G.layout.(f),1);
                    if isempty(jf), error('%s activedim must be either at ''%s'' or ''%s'' location!',f,f,ff), end
                    switch f
                        case 'x'
                            f_off = st.xyoffsets(k,1) + st.xoffset(jf) + sum(st.xoffset(jf+1:end)+st.xstep(jf+1:end));
                            f_step = st.xstep(jf);
                        case 'y'
                            f_off = st.xyoffsets(k,1) + st.yoffset(jf) + sum(st.yoffset(jf+1:end)+st.ystep(jf+1:end));
                            f_step = st.ystep(jf);
                    end
                    
                    % target space between ticks
                    minimumspacing = G.ticksMinimumSpacing(k);
                    fspan = fn_cast(k,st.xspan,st.yspan);
                    minimumspacing = minimumspacing/min(1/fspan(jf),2); % let this target increase up to a factor of two when dimension occupies only a fraction of the space
                    
                    % target space in data coordinates
                    minimumstep = minimumspacing / abs(f_step);
                end
                
                % different display depending on whether header is
                % measure or categorical
                if domeasure
                    [start, scale] = deal(head.start,head.scale);
                    [start, stop] = deal(start, start+(n-1)*scale);
                    if ~dogrid %minimumstep > 2 || ismember(d,G.D.internal_dim)
                        [ticksdata, ticklabels] = G.nicevalues(start,stop,minimumstep*scale);
                        ticksidx = 1 + (ticksdata-start)/scale; % data indices coordinates
                    else
                        % ticks follow data points if minimum step allows
                        % it and we are in an external dimension
                        ticksidx = 1:n;
                        ticksdata = start + (0:n-1)*scale;
                        ticklabels = fn_num2str(ticksdata,'cell');
                    end
                else
                    % ticks for each data point (display only some of
                    % them if there is not enough space for all)
                    ticklabels = row(head.getItemNames());
                    if minimumstep <= 1                     
                        ticksidx = 1:n;
                    else
                        step = G.nicestep(minimumstep);
                        if strcmp(ticklabels{1}, '1')
                            % it seems that we have a mere enumeration,
                            % use a smart step
                            ticksidx = step:step:n;
                        else
                            % make both the first and last appear
                            ticksidx = 1:step:n;
                            if ticksidx(end) ~= n
                                ticksidx(end) = n;
                            end
                        end
                        ticklabels = ticklabels(ticksidx);
                    end
                    % 2D grid: put text rather than using axes ticks
                    % system
                    if dogrid
                        if isempty(st.yspan)
                            row_height = st.yavailable;
                        else
                            row_height = st.yspan(end);
                        end
                        xy = fn_add([0; row_height/2],st.xyoffsets);
                        G.xyticks = gobjects(1,n);
                        for i=1:n
                            G.xyticks(i) = text(xy(1,i),xy(2,i),ticklabels{i}, ...
                                'parent',G.ha,'hittest','off', ...
                                'horizontalalignment','center','verticalalignment','bottom');
                        end
                        ticksidx = []; ticklabels = {};
                    end
                end
                tick = f_off + ticksidx*f_step;
                
                % set ticks!
                if strcmp(f,'x')
                    set(G.ha,'xtick',tick,'xticklabel',ticklabels)
                else
                    set(G.ha,'ytick',fliplr(tick),'yticklabel',fliplr(ticklabels))
                end
            end
            
            % show separations
            if G.showseparation
                G.draw_separations()
            end
            
        end
        function setValueTicks(G)
            % do we show value ticks?
            if ~strcmp(G.D.displaymode,'time courses') ...
                    || (~isempty(G.D.activedim.y) && ~isequal(G.D.activedim.y,G.steps.xydim)) ...
                    || G.D.nodisplay
                ylabel(G.D.ha,'')
                return
            end
            
            % clip values available?
            if isempty(G.D.gridclip)
                error 'programming: setValueTicks should be called after viewdisplay.updateDisplay, so gridclip should be set'
            end
            sz = size(G.D.gridclip); sz(1) = [];
            
            % enough space on y-axis to show values?
            org = G.D.layout;
            st = G.steps;
            minimumspacing = G.ticksMinimumSpacing(2);
            if ~isempty([org.y st.xydim]), minimumspacing = minimumspacing/2; end
%             if st.yavailable < targetspacing
%                 set(G.D.ha,'ytick',[])
%                 ylabel(G.D.ha,'')
%                 return
%             end
            
            % replace min/max clip definitions by min/extent
            startextent = G.D.gridclip;
            startextent(2,:) = startextent(2,:) - startextent(1,:);
            
            % check every 'external' dimension for which we might have
            % several displays on the same y lines (i.e. dimensions at
            % location x, mergeddata or xy/yx); do these display have the same
            % clip values?
            subs0 = substruct('()',repmat({':'},1,1+G.D.nd));
            samecliprow = true;     % time courses on the same row have the same clipping range
            for d = 1:G.D.nd
                if any(d == [org.mergeddata org.x(2:end) G.steps.xydim])
                    % check whether clip values are the same along this
                    % dimension
                    test = diff(startextent,1,1+d);
                    samecliprow = samecliprow && ~any(row(test));
                    if ~samecliprow
                        sameextentrow = ~any(row(test(2,:)));
                        if ~sameextentrow
                            set(G.D.ha,'ytick',[])
                            ylabel(G.D.ha,'')
                        end
                    end
                    % remove this dimension
                    subs = subs0; subs.subs{1+d} = 1;
                    startextent = subsref(startextent,subs);
                elseif any(d == org.y)
                    % nothing to do, keep this dimension
                else
                    % there should be only 1 value in this dimension,
                    % nothing to do
                    if size(startextent,1+d) > 1
                        error 'programming: dimension is not present in the ''external'' layout but gridclip indicates multiple values!'
                    end
                end
            end
            
%             % go back to clip definitions
%             if samecliprow
%                 gridclip(2,:) = gridclip(1,:) + gridclip(2,:);
%             else
%                 gridclip(:,:) = fn_mult([-.5; .5], gridclip(2,:));
%             end
            
            % now gridclip is nonsingleon only in the org.y dimensions;
            % permute dimensions
            startextent = permute(startextent,[1 1+org.y 1+setdiff(1:G.D.nd,org.y)]);
            
            % if not same clip but some extent, clip becomes +/- extent/2
            if ~samecliprow
                startextent(1,:) = -startextent(2,:)/2;
            end
            
            % build ytick and ytickvalues
            if st.xydim
                [nxycol, nxyrow] = deal(st.xyncol,st.xynrow);
            else
                [nxycol, nxyrow] = deal(1);
            end
            ny = prod(sz(org.y));
            [ytick, ytickvalues, yticklabels] = deal(cell([ny nxyrow]));
            for krow = 1:nxyrow
                for ky = 1:ny
                    % middle of the ticks
                    yidx = row(fn_indices(sz(org.y),ky,'g2i'));
                    if ~isempty(org.xy)
                        yrowoffset = st.xyoffsets(2,1+(krow-1)*nxycol);
                    elseif ~isempty(org.yx)
                        yrowoffset = st.xyoffsets(2,krow);
                    else
                        yrowoffset = 0;
                    end
                    yoffset = yrowoffset + sum(st.yoffset) + sum(st.ystep .* yidx);
                    % tick values
                    if any(isnan(startextent(:,ky))), continue, end % happens when data itself consists only of NaNs
                    valuestart = startextent(1,ky);
                    valueextent = startextent(2,ky);
                    minimumsubstep = minimumspacing * valueextent/st.yavailable;
                    [tickvalues, ticklabels] = G.nicevalues(valuestart, valuestart+valueextent, minimumsubstep);
                    yscale = st.yavailable / valueextent;
                    valuemiddle = valuestart + valueextent/2;
                    ytickvalues{ky,krow} = tickvalues;
                    ytick{ky,krow} = yoffset + (tickvalues-valuemiddle)*yscale;
                    yticklabels{ky,krow} = ticklabels;
                end
            end
            % yoffset are descending, so read ytick in reverse order to
            % have only increasing values
            ytick = [ytick{end:-1:1}];
            ytickvalues = [ytickvalues{end:-1:1}];
            yticklabel = [yticklabels{end:-1:1}];
            if ~samecliprow
                % indicate that values are relative to mean
                [yticklabel{ytickvalues==0}] = deal('(mean)');
                idxp = (ytickvalues>0 & ~fn_isemptyc(yticklabel));
                yticklabel(idxp) = fn_map(@(str)['+' str],yticklabel(idxp),'cell');
            end
            
            % set ticks
            set(G.D.ha,'ytick',ytick,'yticklabel',yticklabel)
            ylabel(G.D.ha,'values')
                
            
        end
    end

    % Separations mark for 'external' dimensions
    methods
        function set.showseparation(G, value)
            G.showseparation = value;
            G.draw_separations()
        end
        function set.separationcolor(G, value)
            G.separationcolor = value;
            set(G.separationlines,'color',value)
        end
        function draw_separations(G)
            % no 'smart update', always delete all existing separation
            % lines and redisplay new ones
            deleteValid(G.separationlines)
            G.separationlines = [];
            if ~G.showseparation, return, end
            
            % some properties
            org = G.D.layout;
            nd = G.D.zslice.nd;
            sz = G.D.zslice.sz;
            st = G.steps;
            
            % some variables
            lines = {};
            lwmax = 0; % maximal line-width encountered so far
            
            % x
            kstart = 2;
            linewidth = 0;
            for k = kstart:length(org.x)
                d = org.x(k);                      % dimension for which we will draw vertical lines
                dout = [org.x(k+1:end) st.xydim];  % other dimensions more external than di in x location
                n = sz(d);
                if n == 1, continue, end
                szout = ones(1,nd); szout(dout) = sz(dout);
                nout = prod(szout);
                xpos = zeros(n-1,nout);
                for kout = 1:nout
                    % indices of external dimensions
                    ijk = repmat(fn_indices(szout,kout,'g2i'),[1 n-1]);
                    % indices for dimension d
                    ijk(d,:) = 1.5:n-.5;
                    % x-positions of n-1 lines
                    xpos(:,kout) = sum(fn_add(st.xoffset(k:end)', fn_mult(ijk(org.x(k:end),:),st.xstep(k:end)')),1);
                end
                linewidth = linewidth+.5;
                lwmax = max(lwmax,linewidth);
                lines{end+1} = fn_lines('x',xpos(:),G.D.ha, ...
                    'linewidth',linewidth,'color',G.separationcolor); %#ok<AGROW>
            end
            
            % y
            kstart = 1+strcmp(G.D.displaymode,'image');
            linewidth = 0;
            for k = kstart:length(org.y)
                d = org.y(k);                      % dimension for which we will draw vertical lines
                dout = [org.y(k+1:end) st.xydim];  % other dimensions more external than di in x location
                n = sz(d);
                if n == 1, continue, end
                szout = ones(1,nd); szout(dout) = sz(dout);
                nout = prod(szout);
                ypos = zeros(n-1,nout);
                for kout = 1:nout
                    % indices of external dimensions
                    ijk = repmat(fn_indices(szout,kout,'g2i'),[1 n-1]);
                    % indices for dimension d
                    ijk(d,:) = 1.5:n-.5;
                    % y-positions of n-1 lines
                    ypos(:,kout) = sum(fn_add(st.yoffset(k:end)', fn_mult(ijk(org.y(k:end),:),st.ystep(k:end)')),1);
                end
                linewidth = linewidth+.5;
                lwmax = max(lwmax,linewidth);
                lines{end+1} = fn_lines('y',ypos(:),G.D.ha, ...
                    'linewidth',linewidth,'color',G.separationcolor); %#ok<AGROW>
            end
            
            % xy
            if ~isempty(st.xydim)
                [nxycol, nxyrow] = deal(st.xyncol,st.xynrow);
                xpos = -.5 + (1:ncol-1)/ncol;
                ypos = -.5 + (1:nrow-1)/nrow;
                linewidth = lwmax+.5;
                lines = [lines fn_lines(xpos,ypos,G.D.ha, ...
                    'linewidth',linewidth,'color',G.separationcolor)];
            end
            
            % put lines below other graphic elements
            G.separationlines = fn_map(@row,lines,'array'); % single vector of graphic handles
            uistack(G.separationlines,'bottom')
            
        end
    end
    
    % Coordinates conversions
    methods (Access='private')
        function [subdim, ijk0, mode, invertible] = conversionOptions(G,np,varargin)
            % Options (name/value pairs) for slice/graph conversions:
            % - 'invertible'  [default false] if set to true, exterior
            %               coordinates are rounded, this make the
            %               conversion invertible by calling slice2graph
            % - 'subdim'    [default all dims] dimensions in ijk for which
            %               we perform the conversion; other dimensions
            %               will be assigned to the fixed default values in
            %               ijk0
            % - 'ijk0'      [required if subdim is set] default values for
            %               dimensions where no conversion is requested
            % - 'mode'      value is 'point' or 'vector' [default]

            p = inputParser;
            p.addParameter('subdim',[],@isnumeric)
            p.addParameter('ijk0',[],@isnumeric)
            p.addParameter('mode','point',@(s)ismember(s,{'point', 'vector'}))
            p.addParameter('invertible',false,@islogical)
            parse(p,varargin{:})
            s = p.Results;
            [subdim, ijk0, mode, invertible] = ...
                deal(s.subdim, s.ijk0, s.mode, s.invertible);
            if isempty(subdim)
                subdim = 1:G.D.nd;
            elseif ~isempty(ijk0) && size(ijk0,2)==1 && np>1
                ijk0 = repmat(ijk0,[1 np]);
            end

        end
    end
    methods
        function [zoom, bin] = getZoom(G,varargin)
            % function zoom = getZoom(G[,dim][,'value|effective|indices'])
            % function [offset bin] = getZoom(G[,dim],'off&bin')
            %---
            % Get zoom value for specified dimensions. 
            % Different modes affect the returned value:
            % - 'value' [default]   returns the zooming value specified in
            %               the zoom filter
            % - 'effective' zooming value, after taking into account the
            %               extra space due to not completely filling the
            %               available space (to preserve data aspect ratio)
            % - 'indices'   minimal and maximal displayed data indices
            % - 'displaylimit'   minimal and maximal slice coordinates
            %               being displayed (same as 'indices' for time
            %               courses display, otherwise extends by +/- .5)
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
                case {'indices' 'displaylimit'}
                    zoom = zeros(2,length(dim));
                    for i=1:length(dim), zoom(:,i) = zfilters(i).indicesout([1 end]); end
                    if strcmp(mode,'displaylimit') 
                        if strcmp(G.D.displaymode,'time courses') && ~isempty(G.layout.x)
                            externalxdim = (dim~=G.layout.x(1));
                            zoom(:,externalxdim) = fn_add(zoom(:,externalxdim),[-.5; .5]);                            
                        else
                            zoom = fn_add(zoom,[-.5; .5]);
                        end
                    end
                case 'off&bin'
                    offset = zeros(1,length(dim));
                    for i=1:length(dim), offset(i) = zfilters(i).indicesin(1)-1; end
                    zoom = offset;          % first output: offset
                    bin = [zfilters.bin];   % second output: bin
            end
        end
        function xy = zslice2graph(G,ijk,varargin)
            % function xy = zslice2graph(G,ijk[,options...])
            %---
            % Input:
            % - ijk         index coordinates in the zslice data
            % - options (name/value pairs): see xplr.displaygraph.conversionOptions
            %
            % Output:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            %
            % See also xplr.displaygraph.conversionOptions

            org = G.layout;
            st = G.steps;
            
            % Input points
            np = size(ijk,2);
            [subdim, ijk0, mode, invertible] = conversionOptions(G,np,varargin{:});
            if strcmp(mode,'vector')
                error 'case not handled yet'
            elseif invertible
                warning 'zslice2graph conversion is always invertible, no need to use ''invertible'' flag!'
            end
            if length(subdim) < G.D.nd
                % do not consider all dimensions
                % - for ignored dimensions that are "more exterior" than
                %   dimensions in subdim, we must choose some fixed value
                %   (we choose 1 by default)
                % - ignored dimensions that are "more interior" than 
                %   dimensions in subdim can be totally ignored
                if isempty(ijk0), ijk0 = ones(G.D.nd,np); end
                ijk_ = ijk0;
                if size(ijk,1) == length(subdim)
                    ijk_(subdim,:) = ijk;
                elseif size(ijk,1) == G.D.nd
                    ijk_(subdim,:) = ijk(subdim,:);
                else
                    error('expected %i or %i number of dimensions, but entry points have %i', length(subdim), G.D.nd, size(ijk,1))
                end
                ijk = ijk_;
                orgxok = find(ismember(org.x,subdim),1,'first'):length(org.x); % dimensions on x layout to consider; can be empty
                orgyok = find(ismember(org.y,subdim),1,'first'):length(org.y); % dimensions on x layout to consider; can be empty
            else
                orgxok = 1:length(org.x);
                orgyok = 1:length(org.y);
            end
            
            % "exterior" dimensions must be rounded
            doround = true(1, G.D.nd);
            if ~isempty(orgxok), doround(org.x(orgxok(1))) = false; end
            if ~isempty(orgyok), doround(org.y(orgyok(1))) = false; end
            ijk(doround,:) = round(ijk(doround,:));
            
            x = sum(fn_add(st.xoffset(orgxok)', fn_mult(ijk(org.x(orgxok),:),st.xstep(orgxok)')),1);
            y = sum(fn_add(st.yoffset(orgyok)', fn_mult(ijk(org.y(orgyok),:),st.ystep(orgyok)')),1);
            xy = [x; y];
            if ~isempty(st.xydim)
                xyidx = ijk(st.xydim,:);
                inside = (xyidx>0 & xyidx<=size(st.xyoffsets,2));
                xy(:,inside) = xy(:,inside) + st.xyoffsets(:,xyidx(inside)); 
                % points outside of graph
                xy(:,~inside) = NaN;
            end
        end
        function ijk = graph2zslice(G,xy,varargin)
            % function ijk = graph2zslice(G,xy,options...)
            %---
            % Input:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            % - options (name/value pairs): see xplr.displaygraph.conversionOptions
            %
            % Output:
            % - ijk         index coordinates in the zslice data
            %
            % See also xplr.displaygraph.conversionOptions

            % Input points
            if isvector(xy), xy = xy(:); end
            np = size(xy,2);
            ijk = zeros(G.D.nd,np);
            
            % Parse options
            [subdim, ijk0, mode, invertible] = conversionOptions(G,np,varargin{:});
            if length(subdim)<G.D.nd && isempty(ijk0)
                % Define ijk0 using the first point
                ijk0 = graph2zslice(G,xy(:,1),'invertible',true);
                if np==1, ijk = ijk0(subdim); return, end
            end
            
            % If mode is 'vector', we cannot operate in xy/yx dims, and
            % operate at most on one x and one y dims
            org = G.layout;
            if strcmp(mode,'vector')
                ok = ~any(ismember(subdim,org.xy)) ...
                    && ~any(ismember(subdim,org.yx)) ...
                    && sum(ismember(subdim,org.x)) <= 1 ...
                    && sum(ismember(subdim,org.y)) <= 1;
                if ~ok
                    error 'vector conversion not possible in graph2zslice for this set of dimensions'
                end
            end          
                        
            % xy/yx
            st = G.steps;
            sz = G.zslicesz;
            if strcmp(mode,'point') && ~isempty(st.xydim)
                % take advantage on the fact that the grid spans the full
                % axis
                d = st.xydim;
                x = .5+xy(1,:); y = .5-xy(2,:);
                ncol = round(1/st.xysteps(1)); icol = .5 + x*ncol;
                nrow = round(1/abs(st.xysteps(2))); irow = .5 + y*nrow;
                if ismember(d,subdim)
                    if d==G.layout.xy
                        ijk(d,:) = icol + ncol*round(irow-1);
                    else
                        ijk(d,:) = irow + nrow*round(icol-1);
                    end
                else
                    ijk(d,:) = ijk0(d,:);
                end
                xy = xy - st.xyoffsets(:,fn_coerce(round(ijk(d,:)),1,sz(d)));
            end
            
            % x 
            x = xy(1,:);
            xlayout = org.x;
            for ix = length(xlayout):-1:1
                d = xlayout(ix);
                if strcmp(mode,'point')
                    x = x - st.xoffset(ix);
                end
                if ismember(d,subdim)
                    ijk(d,:) = fn_div(x,st.xstep(ix));
                else
                    ijk(d,:) = ijk0(d,:);
                end
                if strcmp(mode,'point')
                    x = x - fn_mult(round(ijk(d,:)),st.xstep(ix));
                end
            end
            
            % y
            y = xy(2,:);
            ylayout = org.y;
            for iy = length(ylayout):-1:1
                d = ylayout(iy);
                if strcmp(mode,'point')
                    y = y - st.yoffset(iy);
                end
                if ismember(d,subdim)
                    ijk(d,:) = fn_div(y,st.ystep(iy));
                else
                    ijk(d,:) = ijk0(d,:);
                end
                if strcmp(mode,'point')
                    y = y - fn_mult(round(ijk(d,:)),st.ystep(iy));
                end
            end
            
            % we want an output that can be invertible by calling
            % zslice2graph, this means that we should not give the
            % conversion "per dimension" but in a global fashion where
            % dimensions "exterior" to the dimensions of interest are
            % rounded
            if invertible
                doround = true(1, length(ijk));
                if length(subdim) < G.D.nd
                    doround(org.x(find(ismember(org.x,subdim),1,'first'))) = false; % do not round the first dimension of interest on the x location
                    doround(org.y(find(ismember(org.y,subdim),1,'first'))) = false; % do not round the first dimension of interest on the x location
                else
                    if ~isempty(xlayout), doround(xlayout(1)) = false; end
                    if ~isempty(ylayout), doround(ylayout(1)) = false; end
                end
                ijk(doround) = round(ijk(doround));
            end
            
            % not all dimensions?
            if length(subdim) < G.D.nd
                ijk = ijk(subdim,:);
            end
        end
        function xy = slice2graph(G,ijk,varargin)
            % function xy = slice2graph(G,ijk[,options...])
            %---
            % Input:
            % - ijk         index coordinates in the zslice data
            % - options (name/value pairs): see xplr.displaygraph.conversionOptions
            %
            % Output:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            %
            % See also xplr.displaygraph.codnversionOptions

            % Input points
            np = size(ijk,2);
            [subdim, ijk0, mode, invertible] = conversionOptions(G,np,varargin{:});
            if strcmp(mode,'vector')
                error 'case not handled yet'
            elseif invertible
                warning 'zslice2graph conversion is always invertible, no need to use ''invertible'' flag!'
            end
            if length(subdim) < G.D.nd
                if isempty(ijk0)
                    ijk_ = zeros(G.D.nd,np);
                else
                    ijk_ = ijk0;
                end
                if size(ijk,1) == length(subdim)
                    ijk_(subdim,:) = ijk;
                elseif size(ijk,1) == G.D.nd
                    ijk_(subdim,:) = ijk(subdim,:);
                else
                    error('expected %i or %i number of dimensions, but entry points have %i', length(subdim), G.D.nd, size(ijk,1))
                end
                ijk = ijk_;
            end
            
            % first convert from slice to zoomed slice
            if strcmp(mode,'vector')
                error 'case not handled yet'
            else
                [idxoffset, bin] = G.getZoom('off&bin');
                zijk = fn_div(fn_subtract(ijk, idxoffset(:))-.5,bin(:)) + .5;
            end
            
            % then convert to graph coordinates
            if length(subdim)<G.D.nd && isempty(ijk0)
                % let zslice2graph estimate zijk0
                xy = zslice2graph(G,zijk,'mode',mode,'subdim',subdim);
            else
                % ijk and zijk values in other dimensions than subdim have
                % already be assigned using ijk0
                xy = zslice2graph(G,zijk,'mode',mode);
            end
        end
        function ijk = graph2slice(G,xy,varargin)
            % function ijk = graph2slice(G,xy,options...)
            %---
            % Input:
            % - xy          coordinates in the graph (between -0.5 and 0.5)
            % - options (name/value pairs): see xplr.displaygraph.conversionOptions
            %
            % Output:
            % - ijk         index coordinates in the slice data
            %
            % See also xplr.displaygraph.conversionOptions

            % coordinates in zoomed slice
            np = size(xy,2);
            zijk = graph2zslice(G,xy,varargin{:});
            
            % convert to before zooming
            [subdim, ~, mode] = G.conversionOptions(np,varargin{:});
            [idxoffset, bin] = G.getZoom('off&bin');
            if strcmp(mode,'vector')
                ijk = fn_mult(zijk, bin(subdim)');
            else
                ijk = fn_add(idxoffset(subdim)', .5+fn_mult(zijk-.5,bin(subdim)'));
            end
        end
	end

	% Specialized position functions
	methods
        function M = gettransform(G,ijk)
            % function M = gettransform(G,ijk,ybase)
            %---
            % Matrix transformation to place curve/image at ijk data
            % coordinates; note that only coordinates not belonging to the
            % curve/image will be taken into account.
            % 
            % For 'image' displaymode, M will transform indices
            % 1:size(im,1) and 1:size(im,2) to accurate pixel positions.
            % For 'time courses' displaymode, M will transform indices
            % 1:length(x) to accurate x-ordinates, and data values 0 and 1
            % respectively to the bottom and top y-ordinates of the
            % available space.
            % 
            % Input:
            % - ijk     nd * npoint array
            
            st = G.steps;
            org = G.layout;

            % Initialize matrix
            ntransform = size(ijk,2);
            M = repmat(eye(4),[1 1 ntransform]);
            
            % Scale & offset
            % (x)
            if isempty(org.x)
                % no data dimension on x-axis: only 1 data point for time
                % courses or image display, which must be positionned in
                % the center of the available space
                xscale = st.xavailable;
                % index 1 must be positionned at x-ordinate 0, i.e. xoffset + 1*xscale = 0
                xoffsets = -xscale * ones(1,ntransform); 
            else
                xscale = st.xstep(1);
                xoffsets = fn_add( sum(st.xoffset), sum(fn_mult(column(st.xstep(2:end)),ijk(org.x(2:end),:)),1) );
            end
            % (y)
            switch G.D.displaymode
                case 'image'
                    % not possible to have negative scaling in the
                    % hgtransform matrix
                    % -> orienting the images downward will be achieved by
                    % inverting y coordinates at the stage of patch
                    % creation, i.e. will be -1:-1:-size(im,2)
                    if isempty(org.y)
                        % no data dimension on y-axis: similar to above
                        yscale = st.yavailable;
                        % index -1 must be positionned at y-ordinate 0, i.e. yoffset + 1*yscale = 0
                        yoffsets = yscale * ones(1,ntransform);
                    else
                        yscale = abs(st.ystep(1));
                        yoffsets = fn_add( sum(st.yoffset), sum(fn_mult(column(st.ystep(2:end)),ijk(org.y(2:end),:)),1) );
                    end
                case 'time courses'
                    % in this case, the transformation will apply not on
                    % data indices in some dimension, but on data values
                    % -> orient these values upward, and transform them
                    % such that [0 1] fills the available space
                    yscale = st.yavailable;
                    if yscale == 1
                        % occupy full vertical space, because there are no
                        % data dimensions in y, xy or yx location -> get a
                        % nicer display by leaving some gaps above and
                        % below
                        yscale = 1/(1+G.ysep/2);
                    end
                    yoffsets = fn_add( sum(st.yoffset), sum(fn_mult(column(st.ystep(1:end)),ijk(org.y(1:end),:)),1) );
                    % value .5 must be positionned at y-ordinates yoffsets, i.e. yoffsets + .5*yscale = 0
                    yoffsets = yoffsets - yscale/2; 
                otherwise
                    error 'invalid display mode'
            end
            % (add offsets for the xy grid)
            xyoffset = [xoffsets; yoffsets];
            if ~isempty(st.xydim)
                xyoffset = xyoffset + st.xyoffsets(:,ijk(st.xydim,:)); 
            end
            % (set transform matrix)
            M(1,1,:) = xscale;
            M(2,2,:) = yscale; 
            M(1:2,4,:) = xyoffset;
        end
        function [bottomleft, siz] = sub_axes_position(G,ijk)
            % function [bottomleft, siz] = sub_axes_position(G,ijk)
            %---
            % Input:
            % - ijk     nd * npoint array
            % 
            % Output:
            % - pos     4 * npoint array
            
            np = size(ijk,2);
            st = G.steps;

            % Position of sub-axes centers
            % (x)
%             if isempty(org.x)
%                 % no data dimension on x-axis: only 1 data point for time
%                 % courses or image display, which must be positionned in
%                 % the center of the available space
%                 xoffsets = zeros(1,np); 
%             else
                xoffsets = fn_add( sum(st.xoffset(2:end)), sum(fn_mult(column(st.xstep(2:end)),ijk(org.x(2:end),:)),1) );
%             end
            % (y)
            switch G.D.displaymode
                case 'image'
                    yoffsets = fn_add( sum(st.yoffset(2:end)), sum(fn_mult(column(st.ystep(2:end)),ijk(org.y(2:end),:)),1) );
                case 'time courses'
                    yoffsets = fn_add( sum(st.yoffset(1:end)), sum(fn_mult(column(st.ystep(1:end)),ijk(org.y(1:end),:)),1) );
                otherwise
                    error 'invalid display mode'
            end
            % (add offsets for the xy grid)
            xyoffsets = [xoffsets; yoffsets];
            if ~isempty(st.xydim)
                xyoffsets = xyoffsets + st.xyoffsets(:,ijk(st.xydim,:)); 
            end
            
            % We are done!
            siz = repmat([st.xspan(1); st.yspan(1)],[1 np]);
            bottomleft = xyoffsets - siz/2;
        end
        function pos = labelPosition(G,dim,orgin)
            % function pos = labelPosition(G,d[,orgin])
            
            dim = G.D.slice.dimensionNumber(dim);
        
            % steps
            if nargin==3
                st = computeStepsPrivate(G,orgin);
            else
                orgin = G.layout;
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
		function [polygon, center] = selectionMark(G,dim,sel)
            % Create the polygon to display corresponding to a given
            % selection. This is a complex function as it handles many
            % different cases whether the selection is 1D or 2D, which
            % dimensions the selection applies to, and where they are
            % located.
            
            org = G.layout;
            
			% checks
			nd = length(dim);
            dim = G.D.slice.dimensionNumber(dim);
			if sel.nd ~= nd, error 'selection has incorrect number of dimensions', end
            
            % default polygon is empty (no display)
            polygon = nan(2,1); 
            center = nan(2,1); % out of display

			switch nd
				case 1
					lines = sel.polygon; % 2*n array: set of lines
					nline = size(lines,2);
                    % remove lines that are completely out of current view
					zoom = G.getZoom(dim,'value');
                    lines(:, lines(1,:)>zoom(2) | lines(2,:)<zoom(1)) = [];
                    if isempty(lines), return, end
                    % lines spanning beyond the left or right side
                    beyondleft = lines(1,:) < zoom(1);
                    beyondright = lines(2,:) > zoom(2);
                    % clip lines to current view
                    lines(1,beyondleft) = zoom(1);
                    lines(2,beyondright) = zoom(2);
                    % convert from slice to zslice coordinates
                    [idxoffset, bin] = G.getZoom(dim, 'off&bin');
                    lines = (lines-idxoffset-.5)/bin + .5;
                    % display selections as rectangles (for 'x' and 'y'
                    % locations), or as more complex polygon (for 'xy' and
                    % 'yx')
                    st = G.steps;
                    dim_location = G.D.layoutID.dim_locations{dim};
					if ismember(dim_location, {'x' 'y'})
						% convert from zslice to graph coordinates:
						% ignore dimensions that are more internal than dim
						% take value 1 for dimensiont that are more external than dim
						dim_layout = org.(dim_location);
						idx_dim = find(dim_layout==dim,1);
						switch dim_location
							case 'x'
								lines = sum(st.xoffset(idx_dim:end)) + lines*st.xstep(idx_dim) + sum(st.xstep(idx_dim+1:end));
							case 'y'
								lines = sum(st.yoffset(idx_dim:end)) + lines*st.ystep(idx_dim) + sum(st.ystep(idx_dim+1:end));
						end
                        if ~isempty(st.xydim)
                            graphdim = fn_switch(dim_location,'x',1,'y',2);
                            lines = lines + st.xyoffsets(graphdim,1); 
                        end
                        % construct polygon as union of rectangles
                        polygon = cell(1,2*nline-1);
                        for i = 1:nline
                            switch 2*beyondleft(i) + beyondright(i)
                                case 0
                                    % segment within view: full rectangle
                                    polygon{2*i-1} = [lines([1 2 2 1 1],i)'; -.5 -.5 .5 .5 -.5];
                                case 1
                                    % 'rectangle' open on the right side
                                    polygon{2*i-1} = [lines([2 1 1 2],i)'; -.5 -.5 .5 .5];
                                case 2
                                    % 'rectangle' open on the left side
                                    polygon{2*i-1} = [lines([1 2 2 1],i)'; -.5 -.5 .5 .5];
                                case 3
                                    % 'rectangle' open on both sides: 2
                                    % lines
                                    polygon{2*i-1} = [lines([1 2],:)' NaN lines([2 1],i)'; ...
                                        -.5 -.5 NaN .5 .5];
                            end
                        end
                        [polygon{2:2:end}] = deal([NaN; NaN]);
                        polygon = [polygon{:}];
                        center = [mean(lines(:)); 0];
                        
                        % invert coordinates if dim location is 'y'
                        if strcmp(dim_location,'y')
                            polygon = polygon([2 1],:);
                            center = center([2 1]);
                        end
                    elseif strcmp(dim_location, 'xy')
                        % fancy display of selections in grid !
                        hh = abs(st.xysteps(2)) *.46; % half height of the frame around a single grid element
                        polygon = cell(1,2*nline-1);
                        for i = 1:nline
                            % corners: convert from zslice to graph coordinates
                            rline = round(lines(:,i)+[1; -1]*.01);
                            c = st.xyoffsets(:,rline); % 2*2, i.e. x/y * start/stop
                            c(1,:) = c(1,:) + st.xysteps(1) * (lines(:,i) - rline)';
                            % sub-polygon
                            singlerow = (diff(c(2,:)) == 0);
                            if singlerow
                                % the easy case: a simple rectangle
                                % spanning a single line in the grid
                                polygon{2*i-1} = [c(1,[1 2 2 1 1]); c(2,[1 1 2 2 1])+[-1 -1 1 1 -1]*hh];
                            elseif diff(rline) < st.xyncol
                                % two non-intersecting rectangles on two successive lines
                                polygon{2*i-1} = ...
                                    [c(1,1) .5*ones(1,2) c(1,[1 1]) NaN -.5 c(1,[2 2]) -.5*ones(1,2); ...
                                    (c(2,[1 1])+hh) (c(2,[1 1])-hh) (c(2,1)+hh) NaN (c(2,[2 2])+hh) (c(2,[2 2])-hh) (c(2,2)+hh)];
                            else
                                % more difficult: spanning multiple lines
                                hhc = abs(st.xysteps(2)) - hh;
                                polygon{2*i-1} = ...
                                    [c(1,1) .5*ones(1,3) c(1,2)*ones(1,2) -.5*ones(1,3) c(1,1)*ones(1,2); ...
                                    (c(2,[1 1])+hh) (c(2,1)-hhc) (c(2,[2 2])+hhc) (c(2,[2 2])-hh) (c(2,2)+hhc) (c(2,[1 1])-hhc) c(2,1)+hh];
                            end
                        end
                        [polygon{2:2:end}] = deal([NaN; NaN]);
                        polygon = [polygon{:}];
                        
                        % center: will be better pKosition if we average
                        % after conversion from indices to graph positions
                        center = mean(G.slice2graph(sel.dataind,'subdim',dim),2);
                    else
                        error 'not implemented yet'
                        center = [nmean(polygon(1,:)) nmean(polygon(2,:))];
                    end
                case 2
                    % somehow simpler because only the x,y configuration is
                    % allowed
                    % get polygon in slice coordinates
                    polyslice = sel.polygon;
                    centerslice = [mean(polyslice(1,1:end-1)); mean(polyslice(2,1:end-1))]; % remove the last point as it repeats the first one
                    % restrict to the part that is visible within the
                    % current zoom
                    polyslice = G.visiblePolygon(polyslice,dim);
                    display_limit = G.getZoom(dim,'displaylimit');
                    if any(centerslice<display_limit(1,:)' | centerslice>display_limit(2,:)')
                        centerslice(:) = NaN;
                    end
                    % convert to graph
                    polygon = G.slice2graph(polyslice,'subdim',dim);
                    center = G.slice2graph(centerslice,'subdim',dim);
				otherwise
					error 'case not handled yet'
			end

        end
        function output = visiblePolygon(G,polyslice,dim)
                % function output = visiblePolygon(G, polyslice, dim);
                %---
                % @param polyslice: nd*np list of points
                % @param dim: 1*nd list of dimensions in which the
                % polygon is defined
                % @return output: nd*np' list of points: part of the
                % polygon that is within the zoom limits
                % - points that are outside the zoom limits will have their
                % values replaced by NaNs
                % - intermediary points will be inserted exactly at the
                % limit

                % sizes
                [nd, np] = size(polyslice);
                assert(length(dim)==nd)

                % get zoom limits in these dimensions
                zoomSliceValues = G.getZoom(dim,'displaylimit');

                % set the ouput to zeros (they will be set to one if one of the
                % dimension if it's out of display)
                polygonIsOutOfDisplay = false(1,np);
                for dimension = 1:nd  
                   % is equal to one if is out of limits of the zoom or if the
                   % previous value was already 1
                    polygonIsOutOfDisplay = polyslice(dimension,:)<zoomSliceValues(1,dimension) | polyslice(dimension,:)>zoomSliceValues(2,dimension) | polygonIsOutOfDisplay;
                end

                % "Boundary points" will be inserted where there are
                % connections between a displayed point and a point out of
                % display.
                % Example:
                % polygonIsOutOfDisplay = [0  1  1  0  1  0  0  1]
                % boundaryDir =            [1  0 -1  1 -1  0  1]
                % boundaryPrev =            1,    3, 4, 5,    7
                % boundaryNext =            2,    4, 5, 6,    8
                % insertionIndexes =       [2     5  7  8    10] 
                boundaryDir = diff(polygonIsOutOfDisplay);
                boundariesIndexes = find(boundaryDir~=0);
                nboundaries = size(boundariesIndexes,2);
                boundariesValues = nan(nd,nboundaries);
                insertionIndexes = zeros(1,nboundaries);

                % calculate intermediate points between points displayed
                % and points not displayed
                numberOfPointsAdded = 0;
                for k = 1:nboundaries
                    boundaryPrev = boundariesIndexes(k);
                    boundaryNext = boundaryPrev + 1;

                    % which of the two points is inside / outside
                    if(boundaryDir(boundaryPrev)==1)
                        [inside_point, outside_point] = deal(boundaryPrev, boundaryNext);
                    else
                        [inside_point, outside_point] = deal(boundaryNext, boundaryPrev);
                    end

                    % scan dimensions to determine what portion of the
                    % initial segment should be hidden.
                    biggestRatio = 0;
                    for dimension = 1:nd
                        outside_to_inside = polyslice(dimension,inside_point) - polyslice(dimension,outside_point);
                        outside_to_limit_min = zoomSliceValues(1, dimension) - polyslice(dimension,outside_point);
                        outside_to_limit_max = zoomSliceValues(2, dimension) - polyslice(dimension,outside_point);

                        % portion of initial segment must be hidden only if
                        % both values have the same sign (if they have
                        % opposite sign, this means that 'outside_point' is
                        % inside the zoom limit for this dimension)
                        V = [outside_to_limit_max, outside_to_limit_min];
                        if ~any(diff(sign(V)))
                            biggestRatio = max(biggestRatio,min(abs(outside_to_limit_min),abs(outside_to_limit_max))/abs(outside_to_inside));
                        end
                    end

                    % use this ratio to define boundary point
                    boundariesValues(:,k) = polyslice(:,outside_point)+(polyslice(:,inside_point)-polyslice(:,outside_point))*biggestRatio;

                    % insertion position
                    insertionIndexes(k) =  boundaryNext + numberOfPointsAdded;
                    numberOfPointsAdded = numberOfPointsAdded + 1;
                end

                % replace values by NaN for points that must not be
                % displayed
                polyslice(:,polygonIsOutOfDisplay) = NaN;

                % insert boundary points
                output = ones(nd, np + size(boundariesValues,2));
                output(:,setdiff(1:end,insertionIndexes)) = polyslice;
                output(:,insertionIndexes) = boundariesValues;
        end
        function selslice = selection2slice(G,dim,selax)
            % More or less the symmetric of the previous function: convert
            % selection definition in graph coordinates to slice
            % coordinates in the dimensions dim
            % selectionnd object -> use appropriate method for affinity.
            % Works currently only with 2D selections.
            
            if length(dim)~=2, error 'number of dimensions must be 2', end
            
            % use the first point as the origin, work in the zslice to
            % avoid difficulties due to binning
            xy0 = selax.shapes(1).points(:,1);
            zijk0 = round(G.graph2zslice(xy0));

            % infer the affinity matrix graph->zslice
            % (first the linear part)
            xytest = fn_add(xy0, [0 1 0; 0 0 1]);
            zijktest = G.graph2zslice(xytest, 'subdim', dim, 'ijk0', zijk0);
            linearpart = [zijktest(:,2)-zijktest(:,1) zijktest(:,3)-zijktest(:,1)];
            % (then the offset)
            offset = zijktest(:,1) - linearpart * xy0;
            
            % now the affinity matrix graph->slice
            [idxoffset, bin] = G.getZoom('off&bin');
            linearpart = diag(bin(dim)) * linearpart;
            offset = idxoffset(dim)' + .5 + (offset - .5) .* bin(dim)';
            
            % construct affinitynd object
            affinity = xplr.affinitynd(linearpart,offset);
            
            % use the selectionnd method
            selslice = selax.applyaffinity(affinity,G.D.slice.sz(dim));
        end
        %         function [dim, ij] = pointed_dim(G, xy)
        %             % function [dim, ij] = pointed_dim(G, xy)
        %             %---
        %             % Some areas of the main display graph "belong" to different
        %             % dimensions; for example inside an image belongs to the
        %             % dimensions of the images, but outside the image belongs to
        %             % some external dimensions.
        %             % This function is used in particular for zooming with the
        %             % scroll wheel as it helps determining in which dimension to
        %             % perform the zoom. It also returns the coordinated in
        %             % this/these relevant dimensions of the pointed point.
        %             % It returns values for point inside the graph, below it (only
        %             % dimensions in layout.x/.xy/.yx) or on its left (only in
        %             % layout.y/.xy/.yx). It returns an empty array for point
        %             % anywhere else outside the graph.
        %             %
        %             % Input:
        %             % - xy      coordinates of a single point
        %             %
        %             % Output:
        %             % - dim     vector of one or two relevant dimensions for this
        %             %           point
        %             % - ij      coordinate(s) of this point in this(these)
        %             %           dimension(s)
        %
        %             % Check input
        %             assert(isnumeric(xy) && length(xy)==2)
        %
        %             org = G.layout;
        %             st = G.steps;
        %             dim = [];
        %             ij = [];
        %
        %             % Areas not covered (right and top of graph, and the
        %             % bottom-left corner): no dimension selected
        %             if any(xy>.5) || all(xy<-.5)
        %                 return
        %             end
        %
        %         end
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
