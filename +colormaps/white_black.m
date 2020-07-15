function map = white_black(n)
% function map = white_black(n)
%---
% creates a white-black colormap
% 
% See also white_green, white_red

if nargin<1, n=256; end

map = repmat(linspace(1,0,n)',[1 3]);
