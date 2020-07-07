function y = means(varargin)
% function y = means(x1,x2,x3...)
%---
% average all arguments (they need to be same size)

% Thomas Deneux
% Copyright 2004-2017

if nargin<1, help brick.means, return, end

n = 0;
y = [];
for i=1:nargin
    if ~isempty(varargin{i})
        n = n+1;
        xi = varargin{i};
        xi = brick.float(xi);
        if isempty(y)
            y = xi;
        else
            y = y+xi;
        end
    end
end
if n~=1
    y = y/n; 
end