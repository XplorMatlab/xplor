function pr = getPr(var)
%GETPR Get the address of the data stored in a variable: usefull to understand when two variables share the same data in memory
%---
% function pr = getPr(var)
%---
% This MEX-file returns the pointer to the data of variable 'var'. 
% This can be very useful to understand when two Matlab variables having
% the same content do share the same pointer or not, i.e. when is the data
% of these variables shared or duplicated. Ensuring that data of large
% arrays is shared rather than duplicated will ensure efficient memory
% management.

% Thomas Deneux
% Copyright 2015-2017

error 'No MEX-file for brick.getPr was found for your system. Please compile brick.getPr.cpp'
