function map = white_green(n)
% function map = white_green(n)
%---
% creates a white-green colormap
% 
% See also green, black_green

if nargin<1, n=256; end

map = [linspace(1,0,n)' linspace(1,.5,n)' linspace(1,0,n)'];
