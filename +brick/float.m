function x = float(x, flag)
%FLOAT Convert integer to single, keep single or double as such
%---
% function x = float(x[, 'image'])
%---
% convert integer x to single-precision floating number, but do not change
% the class of floating number (in particular double-precision floating
% number remain the same)
% I 'image' flag is used and input is integer, divides values to rescale
% between 0 and 1.

% Thomas Deneux
% Copyright 2012-2017

if nargin==0, help brick.float, return, end

if ~(isnumeric(x) || islogical(x))
    error 'input must be integer'
end
input_type = class(x);
switch input_type
    case {'single' 'double' 'exch.ndSparse' 'ndSparse'}
    case {'int64' 'uint64'}
        x = double(x);
    otherwise
        x = single(x);
end
if nargin>=2
    assert(strcmp(flag, 'image'))
    switch input_type
        case {'single' 'double' 'exch.ndSparse' 'ndSparse'}
        case 'uint8'
            x = x / 255;
        case 'uint16'
            x = x / 65536;
        otherwise
            error 'not handled'
    end
end
