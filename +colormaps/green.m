function map = green(n)
% function map = green(n)
%---
% creates a black-green-white colormap
%
% See also white_green, black_green

if nargin<1, n=256; end
p = floor(n/2);
q = n-p;

blackgreen = [zeros(p,1) (0:p-1)'/(p-1) zeros(p,1)];
greenwhite = [(0:q-1)'/(q-1) ones(q,1) (0:q-1)'/(q-1)];

map = [blackgreen ; greenwhite];