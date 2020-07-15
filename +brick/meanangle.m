function phi = meanangle(angles)
%MEANANGLE Average of angles (result in [-pi pi])
%---
% function phi = meanangle(angles)

% Thomas Deneux
% Copyright 2011-2017

if ~isvector(angles), error('''angles'' must be a vector'), end
z = sum(exp(1i*angles));
phi = angle(z);
