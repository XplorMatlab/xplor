function map = blue(n)
% function map = blue(n)
%---
% creates a black-blue-white colormap

if nargin<1, n=256; end
p = floor(n/2);
q = n-p;

blackblue = [zeros(p,1) zeros(p,1) (0:p-1)'/(p-1)];
bluewhite = [(0:q-1)'/(q-1) (0:q-1)'/(q-1) ones(q,1)];

map = [blackblue ; bluewhite];