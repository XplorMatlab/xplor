function y=subtract(u,v)
% function y=subtract(u,v)
%----
% tool to subtract a matrix row- or column-wise
% ex: y = brick.subtract(rand(3,4),(1:3)')
%     y = brick.subtract(rand(5,2,5),shiftdim(ones(5,1),-2))
%
% See also brick.add, brick.mult, brick.div

% Thomas Deneux
% Copyright 2012-2017

y = bsxfun(@minus,u,v);


    