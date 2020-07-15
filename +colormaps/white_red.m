function map = white_red(n)
% function map = white_red(n)
%---
% creates a white-red colormap
% 
% See also red, black_red

if nargin<1, n=256; end

map = [ones(n,1) (n-1:-1:0)'/(n-1) (n-1:-1:0)'/(n-1)];
