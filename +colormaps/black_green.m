function map = black_green(n)
% function map = black_green(n)
%---
% creates a black-green colormap
% 
% See also green, black_green

if nargin<1, n=256; end

map = [zeros(n,1) linspace(0,1,n)' zeros(n,1)];
