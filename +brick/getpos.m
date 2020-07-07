function pos = getpos(hobj,unit)
%GETPOS Get object position in a specific unit
%---
% function pos = getpos(hobj,unit)
%---
% get the position of specified object according to specific unit

% Thomas Deneux
% Copyright 2015-2017

sunit = get(hobj,'units');
set(hobj,'units',unit)
pos = get(hobj,'pos');
set(hobj,'units',sunit)
