function map = white_blue(n)
% function map = white_blue(n)
%---
% creates a white-red colormap
% 
% See also blue, black_blue

if nargin<1, n=256; end

map = [(n-1:-1:0)'/(n-1) (n-1:-1:0)'/(n-1) ones(n,1)];
