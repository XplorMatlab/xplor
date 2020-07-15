function y = maporient(m)
%MAPCLIP Variation of GRAY colormap with green at the bottom and red at the top.
%       MAPCLIP(M) returns an M-by-3 matrix containing the colormap.
%       MAPCLIP, by itself, is the same length as the current colormap.
%
%       See also mapclip, mapcliplow.

if nargin < 1, m = 129; end
map9 = [1 0 0; 1 1 0; 0 1 0; 0 1 .7; 0 1 1; 0 .5 1; 0 0 1; .5 0 .5; 1 0 0];
y = interp1(linspace(0,1,9),map9,linspace(0,1,m));
