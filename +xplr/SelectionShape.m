classdef SelectionShape
    % selectionshape('point1D',[x1 x2 ...])
    % selectionshape('line1D',[a1 b1 a2 b2 ...])
    % selectionshape('point2D|poly2D|openpoly2D',[x1 x2 ...; y1 y2 ...])
    % selectionshape('rect2D',[x y w h])
    % selectionshape('ellipse2D',{[xc yc],r})            (circle) coordinates of center + radius
    % selectionshape('ellipse2D',{[xc yc],[xu yu],e})    (ellipse) center + principal verctor + eccentricity
    % selectionshape('ring2D',{[xc yc],[xu yu],e,r})     + relative radius
    % selectionshape('ring2D',{[xc yc],[xu yu],[e r]})
    % selectionshape('empty1D|all2D|...')
    % selectionshape('indices',{sizes,[idx1 idx2 ...]})
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
        function S = SelectionShape(type, data)
            if nargin == 0
                return
            end
            
            S.type = type;
            
            switch type
                case 'point1D'
                    S.points = data(:)';
                    S.vectors = zeros(1, 0);
                case 'line1D'
                    lines = brick.row(data);
                    if mod(length(lines), 2) || any(diff(lines(:)) <= 0)
                        error('set of lines should come as ordered non-intersecting segments')
                    end
                    n_line = length(lines)/2;
                    lines = reshape(lines, [1, 2, n_line]);
                    S(n_line).points = []; % pre-allocate
                    for i=1:n_line
                        S(i).points = lines(:, :, i);
                        S(i).vectors = zeros(1, 0);
                    end
                case {'point2D', 'openpoly2D'}
                    if size(data,1) == 1; data = data'; end
                    if size(data,1) ~= 2, error('data should have 2 rows'), end
                    S.points = data;
                    S.vectors = zeros(2, 0);
                case 'poly2D'
                    if size(data,1) ~= 2, error('data should have 2 rows'), end
                    [px, py] = deal(data(1,:), data(2,:)); % poly2cw(data(1,:),data(2,:));
                    S.points = [px(:)'; py(:)'];
                    S.vectors = zeros(2, 0);
                case 'rect2D'
                    if ~isvector(data) || length(data) ~= 4, error('data should be a 4-element vector (x,y,w,h)'), end
                    S.points = [data(1); data(2)];
                    S.vectors = [data(3); data(4)];
                case {'ellipse2D', 'ring2D'}
                    if ~iscell(data)
                        error 'ellipse or ring description must be a cell array';
                    end
                    % center and radius
                    if length(data) < 2 || numel(data{1}) ~= 2 || numel(data{2}) > 2
                        error 'center or radius is ill-defined'
                    end
                    [c, u] = deal(data{1:2});
                    c = c(:);
                    u = u(:);
                    if isscalar(u), u = [u; 0]; end
                    % eccentricity and ring secondary radius
                    switch length(data)
                        case 2
                            if strcmp(type, 'ring2D'), error 'eccentricity and secondary radius missing for ring description', end
                            e = 1;
                            logic = e;
                        case 3
                            logic = data{3};
                        case 4
                            [e, r] = deal(data{3:4});
                            logic = [e, r];
                        otherwise
                            error 'too many elements in shape description'
                    end
                    if length(logic) ~= brick.switch_case(type, 'ellipse2D', 1, 'ring2D', 2)
                        error 'wrong shape definition'
                    end
                    % that's it
                    S.points = c;
                    S.vectors = u;
                    S.logic = logic;
                    S.special = ellipse_vector_to_sym(u, logic(1));
                case 'all'
                    % nothing more needs to be set!
                case 'indices'
                    [sz, indices] = deal(data{:});
                    indices(indices<=0 | indices>prod(sz)) = []; % remove invalid indices
                    S.special = {sz, indices};
                otherwise
                    error('unknown type ''%s''', type)
            end
        end
    end
    
    % Union
    methods
        function S2 = simplify(S)
            % special
            if any(strcmp({S.type}, 'all'))
                S2 = xplr.SelectionShape('all');
                return
            end
            
            % Shapes that cannot be merged
            idx = brick.ismemberstr({S.type}, {'openpoly2D', 'rect2D', 'ellipse2D', 'ring2D'});
            S2 = S(idx);
            S(idx) = [];
            
            % Shapes that are easier to merge at once
            idx = strcmp({S.type}, 'point1D');
            if ~isempty(idx)
                S2 = [S2, xplr.SelectionShape('point1D', unique([S(idx).points]))]; 
            end
            S(idx) = [];
            idx = strcmp({S.type}, 'point2D');
            if ~isempty(idx)
                S2 = [S2, xplr.SelectionShape('point2D', unique([S(idx).points]))]; 
            end
            S(idx) = [];
            idx = strcmp({S.type}, 'line1D');
            if ~isempty(idx)
                S2 = [S2, xplr.SelectionShape('line1D', convert_line_1D(S(idx)))];
            end
            S(idx) = [];
            idx = strcmp({S.type}, 'indices');
            if ~isempty(idx)
                data = cat(1, S(idx).special); % 2 x n
                n = size(data, 2);
                if n == 1
                    S2 = [S2, S(idx)];
                else
                    % check size compatibility
                    sz = data{1, 1}; 
                    for i = 2:n
                        if ~isequal(data{1, i}, sz), error 'indices shape to combine do not have the same size', end
                    end
                    % merge
                    S2 = [S2, xplr.SelectionShape('indices', {sz, unique([data{2, :}])})];
                end
            end
            S(idx) = [];
            
            % Shapes that are easier to merge one after each other
            for k = 1:length(S)
                Sk = S(k);
                f = find(strcmp({S2.type}, Sk.type));
                if isempty(f)
                    S2 = [S2, Sk]; %#ok<AGROW>
                    continue
                end
                switch Sk.type
                    case 'poly2D'
                        % This should be improved, as intersecting
                        % polygons will not be merged
                        S(f).points = [S(f).points, NaN(2, 1), poly2.points];
                    otherwise
                        error('unhandled shape type ''%s''', Sk.type)
                end
            end
            
        end
    end
    
    % Compute indices
    methods
        function ind = indices_1D(S, sizes)
            ind = zeros(1, 0, 'uint32');
            for k=1:length(S)
                switch S(k).type
                    case 'point1D'
                        ind = [ind, round(S(k).points)]; %#ok<AGROW>
                    case 'line1D'
                        line = S(k).points;
                        i_start = max(1, ceil(line(1)));
                        i_end   = min(sizes, floor(line(2)));
                        ind = [ind, i_start:i_end]; %#ok<AGROW>
                    case 'all'
                        ind = ':';
                        return
                    case 'indices'
                        ind = [ind, S(k).special];
                end
            end
            ind(ind < 1 | ind>sizes) = [];
            ind = unique(ind);
        end
        function ind = indices_2D(S, sizes)
            if any(strcmp(S.type, 'all'))
                ind = ':';
                return
            end
            sizes = brick.row(sizes);
            mask = false(sizes);
            for k=1:length(S)
                points = S(k).points;
                vects  = S(k).vectors;
                switch S(k).type
                    case 'point2D'
                        np = size(points, 2);
                        ij = round(points);
                        bad = any(ij < 1 | ij > repmat(sizes(:), 1, np));
                        ij(:, bad) = [];
                        mask(ij(1, :) + sizes(1)*(ij(2, :)-1)) = true;
                    case 'openpoly2D'
                        line = points;
                        longueur = [0, cumsum(sqrt(sum(diff(line, 1, 2).^2)))];
                        [longueur, f] = unique(longueur);
                        line = line(:, f);
                        longueur2 = [0:.1:longueur(end)-.05, longueur(end)];
                        line2 = interp1(longueur, line', longueur2, 'pchip')';
                        np = size(line2, 2);
                        ij = round(line2);
                        bad = any(ij < 1 | ij > repmat(sizes(:), 1, np));
                        ij(:, bad) = [];
                        mask(ij(1, :)+sizes(1)*(ij(2, :)-1)) = true;
                    case 'poly2D'
                        pp = my_poly_split(points);
                        for i=1:length(pp)
                            mask = mask | brick.poly2mask(pp{i}, sizes);
                        end
                    case 'rect2D'
                        i_start = max(1, ceil(points(1)));
                        i_end   = min(sizes(1), floor(points(1)+vects(1)));
                        j_start = max(1, min(ceil(points(2)), round(points(2)+vects(2))));
                        j_end   = min(sizes(2), max(floor(points(2)+vects(2)), round(points(2))));
                        mask(i_start:i_end, j_start:j_end) = true;
                    case 'ellipse2D'
                        points = convert_poly_2D(S(k));
                        mask = mask | brick.poly2mask(points, sizes);
                    case 'ring2D'
                        points = convert_poly_2D(S(k));
                        pp = my_poly_split(points);
                        mask = mask | xor(brick.poly2mask(pp{1}, sizes), brick.poly2mask(pp{2}, sizes));
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
        function lines = convert_line_1D(S)
            lines = zeros(2, 0);
            for k=1:length(S)
                switch S(k).type
                    case 'line1D'
                        lines(:, end+1) = S(k).points; %#ok<AGROW>
                    case {'point1D', 'indices'}
                        if strcmp(S(k).type, 'point1D')
                            indices = sort(round(S(k).points));
                        else
                            indices = S(k).special{2};
                        end
                        if isempty(indices), continue, end
                        gaps = diff(indices) > 1;
                        start = indices([true, gaps]);
                        stop = indices([gaps, true]);
                        lines_k = brick.add(double([start; stop]), [-.5; .5]);
                        lines = [lines, lines_k]; %#ok<AGROW>
                    case 'all'
                        disp 'converting ''all1D'' selection by a [-1e30, 1e30] line'
                        lines = [-1; 1]*1e30;
                        return
                    otherwise
                        error programming
                end
            end
            merge = xor( bsxfun(@lt,lines(1,:)', lines(2,:)), bsxfun(@lt, lines(2,:)', lines(1,:)) );
            for i=1:size(lines,2), merge(i,i) = false; end
            while any(merge(:))
                [i, j] = find(merge, 1, 'first');
                lines(:,i) = [min(lines(1,[i j])), max(lines(2,[i j]))];
                merge(i,:) = merge(i,:) | merge(j,:);
                merge(:,i) = merge(:,i) | merge(:,j);
                merge(i,i) = false;
                lines(:,j) = [];
                merge(j,:) = [];
                merge(:,j) = [];
            end
        end
        function poly = convert_poly_2D(S)
            poly = zeros(2,0);
            for k = 1:length(S)
                Sk = S(k);
                switch Sk.type
                    case 'poly2D'
                        % repeat start point for every component!!
                        idxnan = [find(isnan(Sk.points(1,:))) size(Sk.points,2)+1];
                        n_component = length(idxnan);
                        poly_k = NaN(2, size(Sk.points,2)+n_component);
                        offset = 0;
                        for k = 1:n_component
                            nk = idxnan(k) - offset;
                            poly_k(:, offset+k-1+(1:nk)) = Sk.points(:, offset+[(1:nk-1) 1]);
                            offset = offset + nk;
                        end
                    case 'point2D'
                        p = round(Sk.points);
                        np = size(p,2);
                        poly_k = kron(p,ones(1,6)) + repmat([-.5, -.5, .5, .5, -.5, NaN; -.5, .5, .5, -.5, -.5, NaN], 1, np);
                        poly_k(:,end) = []; % remove last NaN
                    case 'rect2D'
                        poly_k = brick.add(Sk.points, brick.mult(Sk.vectors, [0, 0, 1, 1, 0; 0, 1, 1, 0, 0]));
                    case {'ellipse2D', 'ring2D'}
                        c = Sk.points;
                        u = Sk.vectors;
                        e = Sk.logic(1);
                        phi = linspace(0, 2*pi, 50);
                        u_data = cos(phi);
                        v_data = e*sin(phi);
                        if strcmp(Sk.type, 'ring2D')
                            r = Sk.logic(2);
                            u_data = [u_data, NaN, r*u_data];
                            v_data = [v_data, NaN, r*v_data];
                        end
                        poly_k = brick.add(c, brick.mult(u,u_data) + brick.mult([u(2);-u(1)], v_data));
                    case 'openpoly2D'
                        poly_k = Sk.points;
                    case 'indices'
                        [sz, indices] = deal(Sk.special{:});
                        % first make a mask
                        mask = false([sz, 1]);
                        mask(indices) = true;
                        % then convert to polygon
                        poly_k = brick.mask2poly(mask);
                    case 'all'
                        disp 'converting ''all2D'' selection by a [-1e30, 1e30] square'
                        poly_k = [-1, -1, 1, 1, -1; -1, 1, 1, -1, -1]*1e30;
                    otherwise
                        error programming
                end
                if k == 1
                    poly = poly_k;
                else
                    poly = [poly, [NaN; NaN], poly_k];
                end
            end
        end
    end

    % Affinity
    methods
        function S = apply_affinity(S, affinity)
            if ~isa(affinity, 'xplr.AffinityND')
                error('argument ''afinity'' is expected to be an affinitynd instance')
            end
            for k=1:length(S)
                Sk = S(k);
                S(k).points = affinity.move_points(Sk.points);
                switch Sk.type
                    case 'indices'
                        error 'affinity cannot be performed on selection of type ''indices'''
                    case {'ellipse2D', 'ring2D'}
                        % this is a bit bad that ellipse main axis vector
                        % and eccentricity cannot be dealt like usual
                        % vectors and logic
                        [S(k).vectors S(k).logic S(k).special] = ...
                            ellipse_affinity(Sk.vectors, Sk.logic(1), Sk.special, affinity.linear_part);
                        if strcmp(Sk.type, 'ring2D'), S(k).logic(2) = Sk.logic(2); end
                    otherwise
                        S(k).vectors = affinity.move_vectors(Sk.vectors);
                        % 'logic' remains unchanged
                end
            end
        end
    end
    
    % Misc
    methods
        function b = is_point(S, tol)
            switch S.type
                case {'point1D', 'point2D', 'line1D', 'openpoly2D', 'poly2D'}
                    b = all(all(abs(diff(S.points, 1, 2)) <= tol));
                case 'rect2D'
                    b = all(abs(S.vectors) <= tol);
                case {'ellipse2D' 'ring2D'}
                    b = all(abs(S.vectors) <= tol);
                case 'all'
                    b = false;
                case 'indices'
                    b = (length(S.special{2}) <= 1);
                otherwise
                    error programming
            end
        end
        function S2 = to_point(S)
            if strfind(S.type, '1D')
                S2 = xplr.SelectionShape('point1D', mean(S.points));
            elseif strfind(S.type, '2D')
                S2 = xplr.SelectionShape('point2D', mean(S.points,2));
            else
                error 'not handled yet'
            end
        end
    end
end


%------------
% ELLIPSE
%------------

% Ellipse defined either by center, main radius vector and eccentricity, or
% by its bilinear equation (x-c)'A(x-c) = 1,

function [u, e, A] = ellipse_affinity(u, e, A, M)

    % ellipse equation becomes, for y=Mx: (y-Mc)'(M^-1' A M^-1)(y-Mc) = 1
    M1 = M^-1;
    A = M1'*A*M1;
    [u, e] = ellipse_sym_to_vector(A);

end

%---
function A = ellipse_vector_to_sym(u,e)

    % (U,c) is the referential of the ellipse, x->y=U'(x-c) returns coordinates
    % in this referential, in which the ellipse equation is
    % y(1)^2 + y(2)^2/e^2 = r^2
    r = norm(u);
    u = u/r;
    U = [[-u(2); u(1)], u];
    A = U*(diag([1/e^2 1])/r^2)*U';

end

%---
function [u, e] = ellipse_sym_to_vector(A)

    % eigenvalue decomposition returns U, r and e as above
    [U, D] = svd(A); % better use svd than eig because output is real even if A is not exactly symmetric
    r = D(2,2)^-(1/2);
    e = D(1,1)^-(1/2) / r;
    u = r*U(:,2);

end

%---
function pp = my_poly_split(points)

    ksep = find(any(isnan(points), 1));
    n = length(ksep) + 1;
    if n == 1
        pp = {points};
    else
        pp = cell(1, n);
        ok = false(1, n);
        ksep = [0, ksep, size(points,2)+1];
        for i=1:n
            pp{i} = points(:, ksep(i)+1:ksep(i+1)-1);
            ok(i) = ~isempty(pp{i});
        end
        pp = pp(ok);
    end

end
