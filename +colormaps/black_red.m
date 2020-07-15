function map = black_red(n)
% function map = black_red(n)
%---
% creates a black-red colormap
% 
% See also red, white_red

if nargin<1, n=256; end

map = [linspace(0,1,n)' zeros(n,1) zeros(n,1)];
