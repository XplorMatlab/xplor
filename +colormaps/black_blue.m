function map = black_blue(n)
% function map = black_blue(n)
%---
% creates a black-blue colormap
% 
% See also blue, white_blue

if nargin<1, n=256; end

map = [zeros(n,1) zeros(n,1) linspace(0,1,n)'];
