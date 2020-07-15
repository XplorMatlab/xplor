function b = isemptyc(varargin)
%ISEMPTYC Which elements of a cell array are empty
%---
% function b = isemptyc(c)
% function b = isemptyc(x1,x2,x3,...)
%---
% returns an array of logicals of the same size as cell array c indicating
% which elements of c are empty
%
% See also brick.find, brick.map, brick.itemlengths

% Thomas Deneux
% Copyright 2011-2017

if nargin==0, help brick.isemptyc, return, end

if nargin==1
    c = varargin{1};
else
    c = varargin;
end

b = false(size(c));
for k=1:numel(c), b(k) = isempty(c{k}); end
