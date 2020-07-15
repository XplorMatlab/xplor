function n = itemlengths(c)
%ITEMLENGTHS Return lengths of all elements inside a cell array
%---
% function n = itemlengths(c)
%---
% returns an array of same size as cell array c, containing the length of
% each of its elements

% Thomas Deneux
% Copyright 2015-2017

n = zeros(size(c));
for i=1:numel(c), n(i) = length(c{i}); end
