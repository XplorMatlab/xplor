function map = signcheck(n)
% function map = signcheck(n)
%---
% creates a 2-colors colormap (one for each sign, gray in the middle)
%
% See also green, mapclip

if nargin<1, n=255; end
p = floor(n/2);

blue2 = [0 0 .6];
blue1 = [.5 .5 1];
grey = [1 1 1]*.7;
yellow1 = [.9 .9 .8];
yellow2 = [1 1 0];

blue = brick.add(blue1, brick.mult((blue2-blue1),(p-1:-1:0)'/(p-1)));
yellow = brick.add(yellow1, brick.mult((yellow2-yellow1),(0:p-1)'/(p-1)));

if mod(n,2)
    map = [blue ; grey ; yellow];
else
    map = [blue ; yellow];
end

if nargout==0
    colormap(map)
    clear map
end