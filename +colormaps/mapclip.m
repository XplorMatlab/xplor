function y = mapclip(m)
%MAPCLIP Variation of GRAY colormap with green at the bottom and red at the top.
%       MAPCLIP(M) returns an M-by-3 matrix containing the colormap.
%       MAPCLIP, by itself, is the same length as the current colormap.
%
%       See also mapcliplow, mapcliphigh.
%
%       Doron, 17 Nov 1994.
%

if nargin < 1, m = size(get(gcf,'colormap'),1); end
y = (0:m-1)'/max(m-1,1);
y = [y y y];

y(1,:) = [0 0 .5];
y(m,:) = [1 0 0];
