function [m I] = max(a)
% function [m I] = max(a)
%---
% find the global max in an array, and give its coordinates
%
% See also brick.min

% Thomas Deneux
% Copyright 2005-2017

if nargin==0, help brick.min, return, end

[m i] = max(a(:));
i = i-1;
s = size(a);
for k=1:length(s)
    I(k) = mod(i,s(k))+1;
    i = floor(i/s(k));
end
