function map = red(n)
% function map = red(n)
%---
% creates a black-red-white colormap
% 
% See also white_red, black_red

if nargin<1, n=256; end
p = floor(n/2);
q = n-p;

blackred = [(0:p-1)'/(p-1) zeros(p,1) zeros(p,1)];
redwhite = [ones(q,1) (0:q-1)'/(q-1) (0:q-1)'/(q-1)];

map = [blackred ; redwhite];