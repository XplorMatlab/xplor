classdef selectionshape
    % selectionshape('point1D',[x1 x2 ...])
    % selectionshape('line1D',[a1 b1 a2 b2 ...])
    % selectionshape('point2D|line2D|poly2D',[x1 x2 ...; y1 y2 ...])
    % selectionshape('rect2D',[x y w h])
    % selectionshape('ellipse2D',{[xc yc],r})            (circle) coordinates of center + radius
    % selectionshape('ellipse2D',{[xc yc],[xu yu],e})    (ellipse) center + principal verctor + eccentricity
    % selectionshape('ring2D',{[xc yc],[xu yu],e,r})     + relative radius
    % selectionshape('ring2D',{[xc yc],[xu yu],[e r]})
    % selectionshape('empty1D|all2D|...')
    % selectionshape('indices',[idx1 idx2 ...],sizes)
    % selectionshape('product',[sel1 sel2 ...])          sel1, sel2, ... are selectionND themselves
    properties
        type
        points      % moved as points by affine transformations
        vectors     % moved as points by affine transformations
        logic       % unchanged by affine transformations
        special     % specific code depending on type for applying affine transformations
    end
    
    % Constructor
    methods
        function S = selectionshape(type,data)
            if nargin == 0
                return
            end
            
            S.type = type;
            
            switch type
                case 'point1D'
                    S.points = data(:)';
                    S.vectors = zeros(1,0);
                case 'line1D'
                    lines = row(data);
                    if mod(length(lines),2) || any(diff(lines(:))<=0)
                        error('set of lines should come as ordered non-intersecting segments')
                    end
                    nline = length(lines)/2;
                    lines = reshape(lines,[1 2 nline]);
                    S(nline).points = []; % pre-allocate
                    for i=1:nline
                        S(i).points = lines(:,:,i);
                        S(i).vectors = zeros(1,0);
                    end
                case 'point2D'
                    if size(data,1)==1; data = data'; end
                    if size(data,1)~=2, error('data should have 2 rows'), end
                    S.points = data;
                    S.vectors = zeros(2,0);
                case 'line2D'
                    if size(data,1)==1; data = data'; end
                    if size(data,1)~=2, error('data should have 2 rows'), end
                    if mod(size(data,2),2), error('set of lines improperly defined'), end
                    nline = size(data,2)/2;
                    lines = reshape(data,[2 2 nline]);
                    S(nline).points = []; % pre-allocate
                    for i=1:nline
                        S(i).points = lines(:,:,i);
                        S(i).vectors = zeros(2,0);
                    end
                case 'poly2D'
                    if size(data,1)~=2, error('data should have 2 rows'), end
                    if ~all(data(:,1)==data(:,end)), data = data(:,[1:end 1]); end
                    [px, py] = deal(data(1,:),data(2,:)); % poly2cw(data(1,:),data(2,:));
                    S.points = [px(:)'; py(:)'];
                    S.vectors = zeros(2,0);
                case 'rect2D'
                    if ~isvector(data) || length(data)~=4, error('data should be a 4-element vector (x,y,w,h)'), end
                    if any(data(3:4)<0), error('width and height must be >= 0'), end
                    S.points = [data(1); data(2)];
                    S.vectors = [data(3); data(4)];
                case {'ellipse2D' 'ring2D'}
                    if ~iscell(data)
                        error 'ellipse or ring description must be a cell array';
                    end
                    % center and radius
                    if length(data)<2 || numel(data{1})~=2 || numel(data{2})>2
                        error 'center or radius is ill-defined'
                    end
                    [c, u] = deal(data{1:2});
                    c = c(:); u = u(:);
                    if isscalar(u), u = [u; 0]; end
                    % eccentricity and ring secondary radius
                    switch length(data)
                        case 2
                            if strcmp(type,'ring2D'), error 'eccentricity and secondary radius missing for ring description', end
                            e = 1;
                            logic = e;
                        case 3
                            logic = data{3};
                        case 4
                            [e r] = deal(data{3:4});
                            logic = [e r];
                        otherwise
                            error 'too many elements in shape description'
                    end
                    if length(logic)~=fn_switch(type,'ellipse2D',1,'ring2D',2)
                        error 'wrong shape definition'
                    end
                    % that's it
                    S.points = c;
                    S.vectors = u;
                    S.logic = logic;
                    S.special = EllipseVector2Sym(u,logic(1));
                case 'all'
                    % nothing more needs to be set!
                case 'indices'
                    S.special = data;
                otherwise
                    error('unknown type ''%s''',type)
            end
        end
    end
    
    % Union
    methods
        function S2 = simplify(S)
            % special
            if any(strcmp({S.type},'all'))
                S2 = xplr.selectionshape('all');
                return
            end
            
            % Shapes that cannot be merged
            idx = fn_ismemberstr({S.type},{'line2D' 'rect2D' 'ellipse2D' 'ring2D'});
            S2 = S(idx);
            S(idx) = [];
            
            % Shapes that are easier to merge at once
            idx = strcmp({S.type},'point1D');
            if ~isempty(idx)
                S2 = [S2 xplr.selectionshape('point1D',unique([S(idx).points]))]; 
            end
            S(idx) = [];
            idx = strcmp({S.type},'point2D');
            if ~isempty(points2D)
                S2 = [S2 xplr.selectionshape('point2D',unique([S(idx).points]))]; 
            end
            S(idx) = [];
            idx = strcmp({S.type},'line1D');
            if ~isempty(lines1D)
                S2 = [S2 xplr.selectionshape('line1D',ConvertLine1D(S(idx)))];
            end
            S(idx) = [];
            idx = strcmp({S.type},'indices');
            if ~isempty(lines1D)
                S2 = [S2 xplr.selectionshape('indices',unique([S(idx).special]))];
            end
            S(idx) = [];
            
            % Shapes that are easier to merge one after each other
            for k = 1:length(S)
                Sk = S(k);
                f = find(strcmp({S2.type},Sk.type));
                if isempty(f)
                    S2 = [S2 Sk]; %#ok<AGROW>
                    continue
                end
                switch Sk.type
                    case 'poly2D'
                        % This should be improved, as intersecting
                        % polygons will not be merged
                        S(f).points = [S(f).points NaN(2,1) poly2.points];
                    otherwise
                        error('unhandled shape type ''%s''', Sk.type)
                end
            end
            
        end
    end
    
    % Compute indices
    methods
        function ind = indices1D(S,sizes)
            ind = zeros(1,0,'uint32');
            for k=1:length(S)
                switch S(k).type
                    case 'point1D'
                        ind = [ind round(S(k).points)]; %#ok<AGROW>
                    case 'line1D'
                        line = S(k).points;
                        istart = max(1,ceil(line(1)));
                        iend   = min(sizes,floor(line(2)));
                        ind = [ind istart:iend]; %#ok<AGROW>
                    case 'all'
                        ind = ':';
                        return
                    case 'indices'
                        ind = [ind S(k).special];
                end
            end
            ind(ind<1 | ind>sizes) = [];
            ind = unique(ind);
        end
        function ind = indices2D(S,sizes)
            if any(strcmp(S.type,'all'))
                ind = ':';
                return
            end
            mask = false(sizes);
            for k=1:length(S)
                points = S(k).points;
                vects  = S(k).vectors;
                switch S(k).type
                    case 'point2D'
                        np = size(points,2);
                        ij = round(points);
                        bad = any(ij<1 | ij>repmat(sizes(:),1,np));
                        ij(:,bad) = [];
                        mask(ij(1,:)+sizes(1)*(ij(2,:)-1)) = true;
                    case 'line2D'
                        line = points;
                        longueur = [0 cumsum(sqrt(sum(diff(line,1,2).^2)))];
                        [longueur f] = unique(longueur);
                        line = line(:,f);
                        longueur2 = [0:.1:longueur(end)-.05 longueur(end)];
                        line2 = interp1(longueur,line',longueur2,'pchip')';
                        np = size(line2,2);
                        ij = round(line2);
                        bad = any(ij<1 | ij>repmat(sizes(:),1,np));
                        ij(:,bad) = [];
                        mask(ij(1,:)+sizes(1)*(ij(2,:)-1)) = true;
                    case 'poly2D'
                        pp = mypolysplit(points);
                        for i=1:length(pp)
                            mask = mask | fn_poly2mask(pp{i},sizes);
                        end
                    case 'rect2D'
                        istart = max(1,ceil(points(1)));
                        iend   = min(sizes(1),floor(points(1)+vects(1)));
                        jstart = max(1,min(ceil(points(2)),round(points(2)+vects(2))));
                        jend   = min(sizes(2),max(floor(points(2)+vects(2)),round(points(2))));
                        mask(istart:iend,jstart:jend) = true;
                    case 'ellipse2D'
                        tmp = ConvertPoly2D(S(k));
                        points = tmp.points;
                        mask = mask | fn_poly2mask(points,sizes);
                    case 'ring2D'
                        tmp = ConvertPoly2D(S(k));
                        pp = mypolysplit(tmp.points);
                        mask = mask | xor(fn_poly2mask(pp{1},sizes),fn_poly2mask(pp{2},sizes));
                    case 'indices'
                        mask(indices) = true;
                    otherwise
                        error programming
                end
            end
            ind = uint32(find(mask));
        end
    end
    
    % Conversion
    methods
        function lines = ConvertLine1D(S)
            lines = zeros(2, 0);
            for k=1:length(S)
                switch S(k).type
                    case 'line1D'
                        lines(:,end+1) = S(k).points; %#ok<AGROW>
                    case {'point1D' 'indices'}
                        if strcmp(S(k).type,'point1D')
                            indices = sort(round(S(k).points));
                        else
                            indices = S(k).special;
                        end
                        if isempty(indices), continue, end
                        gaps = diff(indices)>1;
                        start = indices([true gaps]);
                        stop = indices([gaps true]);
                        linesk = fn_add(double([start; stop]), [-.5; .5]);
                        lines = [lines linesk]; %#ok<AGROW>
                    case 'all'
                        disp 'converting ''all1D'' selection by a [-1e30 1e30] line'
                        lines = [-1; 1]*1e30;
                        return
                    otherwise
                        error programming
                end
            end
            merge = xor( bsxfun(@lt,lines(1,:)',lines(2,:)), bsxfun(@lt,lines(2,:)',lines(1,:)) );
            for i=1:size(lines,2), merge(i,i) = false; end
            while any(merge(:))
                [i, j] = find(merge,1,'first');
                lines(:,i) = [min(lines(1,[i j])) max(lines(2,[i j]))];
                merge(i,:) = merge(i,:) | merge(j,:);
                merge(:,i) = merge(:,i) | merge(:,j);
                merge(i,i) = false;
                lines(:,j) = [];
                merge(j,:) = [];
                merge(:,j) = [];
            end
        end
        function poly = ConvertPoly2D(S)
            poly = zeros(2,0);
            for k = 1:length(S)
                Sk = S(k);
                switch Sk.type
                    case 'poly2D'
                        polyk = Sk.points;
                    case 'point2D'
                        p = round(Sk.points);
                        np = size(p,1);
                        polyk = kron(p,ones(1,6)) + repmat([-.5 -.5 .5 .5 -.5 NaN; -.5 .5 .5 -.5 -.5 NaN],1,np);
                        polyk(end,:) = []; % remove last NaN
                    case 'rect2D'
                        polyk = fn_add(Sk.points,fn_mult(Sk.vectors,[0 0 1 1 0; 0 1 1 0 0]));
                    case {'ellipse2D' 'ring2D'}
                        c = Sk.points;
                        u = Sk.vectors;
                        e = Sk.logic(1);
                        phi = linspace(0,2*pi,20);
                        udata = cos(phi);
                        vdata = e*sin(phi);
                        if strcmp(Sk.type,'ring2D')
                            r = Sk.logic(2);
                            udata = [udata NaN r*udata];
                            vdata = [vdata NaN r*vdata];
                        end
                        polyk = fn_add(c,fn_mult(u,udata)+fn_mult([u(2);-u(1)],vdata));
                    case 'line2D'
                        polyk = Sk.points;
                    case 'all'
                        disp 'converting ''all2D'' selection by a [-1e30 1e30] square'
                        polyk = [-1 -1 1 1 -1; -1 1 1 -1 -1]*1e30;
                    otherwise
                        error programming
                end
            end
            if k == 1
                poly = polyk;
            else
                poly = [poly [NaN; NaN] polyk];
            end
        end
    end

    % Affinity
    methods
        function S = applyaffinity(S,mat)
            for k=1:length(S)
                Sk = S(k);
                xpoints = [ones(1,size(Sk.points,2)); Sk.points];
                ypoints = mat*xpoints;
                S(k).points = ypoints(2:end,:);
                switch Sk.type
                    case 'indices'
                        error 'affinity cannot be performed on selection of type ''indices'''
                    case {'ellipse2D' 'ring2D'}
                        % this is a bit bad that ellipse main axis vector
                        % and eccentricity cannot be dealt like usual
                        % vectors and logic
                        [S(k).vectors S(k).logic S(k).special] = ...
                            EllipseAffinity(Sk.vectors,Sk.logic(1),Sk.special,mat(2:3,2:3));
                        if strcmp(Sk.type,'ring2D'), S(k).logic(2) = Sk.logic(2); end
                    otherwise
                        S(k).vectors = mat(2:end,2:end)*Sk.vectors;
                        % 'logic' remains unchanged
                end
            end
        end
    end
    
    % Misc
    methods
        function b = ispoint(S,tol)
            switch S.type
                case {'point1D','point2D','line1D','line2D','poly2D'}
                    b = all(all(abs(diff(S.points,1,2))<=tol));
                case 'rect2D'
                    b = all(abs(S.vectors)<=tol);
                case {'ellipse2D' 'ring2D'}
                    b = all(abs(S.vectors)<=tol);
                case 'all'
                    b = false;
                otherwise
                    error programming
            end
        end
    end
end


%------------
% ELLIPSE
%------------

% Ellipse defined either by center, main radius vector and eccentricity, or
% by its bilinear equation (x-c)'A(x-c) = 1,

function [u, e, A] = EllipseAffinity(u,e,A,M)

    % ellipse equation becomes, for y=Mx: (y-Mc)'(M^-1' A M^-1)(y-Mc) = 1
    M1 = M^-1;
    A = M1'*A*M1;
    [u, e] = EllipseSym2Vector(A);

end

%---
function A = EllipseVector2Sym(u,e)

    % (U,c) is the referential of the ellipse, x->y=U'(x-c) returns coordinates
    % in this referential, in which the ellipse equation is
    % y(1)^2 + y(2)^2/e^2 = r^2
    r = norm(u);
    u = u/r;
    U = [[-u(2); u(1)] u];
    A = U*(diag([1/e^2 1])/r^2)*U';

end

%---
function [u e] = EllipseSym2Vector(A)

    % eigenvalue decomposition returns U, r and e as above
    [U D] = svd(A); % better use svd than eig because output is real even if A is not exactly symmetric
    r = D(2,2)^-(1/2);
    e = D(1,1)^-(1/2) / r;
    u = r*U(:,2);

end

%---
function pp = mypolysplit(points)

    ksep = find(any(isnan(points),1));
    n = length(ksep)+1;
    if n==1
        pp = {points};
    else
        pp = cell(1,n); ok = false(1,n);
        ksep = [0 ksep size(points,2)+1];
        for i=1:n
            pp{i} = points(:,ksep(i)+1:ksep(i+1)-1);
            ok(i) = ~isempty(pp{i});
        end
        pp = pp(ok);
    end

end