function y = add(u,v,varargin)
% function y = add(u,v,...)
%---
% tool to add matrices and vectors
% ex: y = brick.add(rand(4,5),(1:4)')
%     y = brick.add(1:5,(1:4)')
%
% See also brick.mult, brick.subtract, brick.div, brick.eq

% Thomas Deneux
% Copyright 2002-2017

y = bsxfun(@plus,u,v);
for i=1:length(varargin), y = bsxfun(@plus,y,varargin{i}); end
    
