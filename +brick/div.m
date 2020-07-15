function y=div(u,v)
% function y=div(u,v)
%----
% tool to divide a matrix row- or column-wise
% ex: y = brick.div(rand(3,4),(1:3)')
%     y = brick.div(rand(5,2,5),shiftdim(ones(5,1),-2))
%
% See also brick.mult, brick.add, brick.subtract

% Thomas Deneux
% Copyright 2012-2017

y = bsxfun(@rdivide,u,v);

    