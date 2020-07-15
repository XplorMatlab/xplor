function y=mult(u,v,varargin)
% function y=mult(u,v,...)
%----
% tool to multiply a matrix row- or column-wise
% ex: y = brick.mult(rand(3,4),(1:3)')
%     y = brick.mult(rand(5,2,5),shiftdim(ones(5,1),-2))
%
% See also brick.add, brick.subtract, brick.div

% Thomas Deneux
% Copyright 2002-2017

% Check sizes for first 2 arguments
su = size(u); sv = size(v); 
n = min(length(su),length(sv));
su = su(1:n); sv = sv(1:n);
if ~all(su==sv | su==1 | sv==1)
    error('Incompatible sizes for broadcast operation: [%s] vs. [%s].', ...
        brick.strcat(size(u),' '), brick.strcat(size(v),' '))
end

% Perform broadcast operation
y = bsxfun(@times,u,v);
for i=1:length(varargin), y = bsxfun(@times,y,varargin{i}); end

    