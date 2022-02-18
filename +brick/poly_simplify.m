function poly = poly_simplify(poly, flag)
% Polygon operations by wrapping the polyshape class introduced in R2017b
% p = poly_simplify(p [, tol])
% {p} = poly_simplify({p} [, tol])
% p = poly_simplify({p}, 'union')
% p = poly_simplify({p}, 'intersect')
%---
% Performs polygon simplification, union or intersection.
% For simplification: tol is the tolerance used in rmslivers function,
% default value is 2e-4. Use tol=0 for simply fixing up duplicate vertices
% and other degeneracies.

% Operation specification
if nargin<2
    flag = 'simplify'; 
    tol = 2e-4;
elseif isnumeric(flag)
    tol = flag;
    flag = 'simplify';
end

% Input polygon(s)
if iscell(poly)
    for k = 1:length(poly)
        x = poly{k};
        poly{k} = polyshape(x(1,:), x(2,:));
    end
    poly = [poly{:}];
else
    if ~strcmp(flag, 'simplify')
        error 'for ''union'' and ''intersect'' operations, input polygons must be a cell array'
    end
    poly = polyshape(poly(1,:), poly(2,:));
end

% Operation
switch flag
    case 'simplify'
        if tol>0, poly = rmslivers(poly, tol); end
    case 'union'
        poly = poly.union();
        poly = rmslivers(poly, 1e-6); % small cleanup
    case 'intersect'
        poly = poly.intersect();
        poly = rmslivers(poly, 1e-6); % small cleanup
    otherwise
        error('unknown polygon operation ''%s''', flag)
end

% repeat start point for every component!!

% Output
if isscalar(poly)
    poly = repeat_start_point(poly.Vertices');
else
    poly = brick.map(@(x)repeat_start_point(x.Vertices'), poly);
end
    

%---
function poly2 = repeat_start_point(poly)

if isempty(poly)
    poly2 = poly;
    return
end
idxnan = [find(isnan(poly(1,:))) size(poly,2)+1];
n_component = length(idxnan);
poly2 = NaN(2, size(poly,2)+n_component);
offset = 0;
for k = 1:n_component
    nk = idxnan(k) - offset;
    poly2(:, offset+k-1+(1:nk)) = poly(:, offset+[(1:nk-1) 1]);
    offset = offset + nk;
end



