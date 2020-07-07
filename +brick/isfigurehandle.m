function [b c] = isfigurehandle(x)
%ISFIGUREHANDLE Is handle a plausible figure handle
%---
% function [isfig isnewfig] = isfigurehandle(x)
%---
% This function detects whether x can be a valid figure handle (i.e. is a
% positive integer), and opens the figure if it does not exist yet.

% Thomas Deneux
% Copyright 2008-2017

if ~isscalar(x)
    b = false(size(x));
    for i=1:numel(x), b(i) = brick.isfigurehandle(x(i)); end
    return
end

b = (ishandle(x) && strcmp(get(x,'type'),'figure'));
c = false;
if ~b && isnumeric(x) && mod(x,1)==0 && x>0
    b = true;
    c = true;
    figure(x)
end

if nargout==0
    clear b
end
