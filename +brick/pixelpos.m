function pos = pixelpos(hobj,varargin)
% function pos = pixelpos(hobj[,'recursive'][,'strict'])
%---
% returns the position in pixels of any object without needing to
% change any units values
%
% In R2014b and later, wraps function getpixelposition
%
% Input:
% - hobj            most graphic objects
% - 'recursive'     return position relative to parent figure rather than
%                   to immediate parent
% - 'strict'        take into account that when DataAspectRatio is
%                   fixed, only a subset of the axes space is used
%
% See also brick.pixelsize, getpixelposition

% Thomas Deneux
% Copyright 2011-2017

if nargin==0, help brick.pixelpos, return, end
[recursive, strict] = brick.flags({'recursive' 'strict'},varargin);
 
% get pixel position
if strcmp(get(hobj,'type'),'text')
    % special: get(hobj,'pos') returns a 3-element vector, and
    % getpixelposition(hobj) would return [0 0 0 0]
    curunit = get(hobj,'unit');
    set(hobj,'unit','pixel')
    pos = get(hobj,'pos'); pos(3:4) = 0;
    set(hobj,'unit',curunit)
    if recursive
        ppos = brick.pixelpos(get(hobj,'parent'),'recursive','strict');
        pos(1:2) = (ppos(1:2)-1) + pos(1:2);
    end
else
    pos = getpixelposition(hobj,recursive);
end

% 'strict' flag: take into account the fact that when DataAspectRatio is
% fixed, only a subset of the space is used
if strict && strcmp(get(hobj,'type'),'axes') && strcmp(get(hobj,'dataaspectratiomode'),'manual')
    pos0 = pos(1:2);
    psiz = pos(3:4);
    % first find out which dimension is not fully occupied
    availableratio = psiz(1)/psiz(2);
    dataratio = get(hobj,'dataaspectratio'); % ratio(1) in x should be the same length as ratio(2) in y
    actualratio = (diff(get(hobj,'xlim'))/dataratio(1)) / (diff(get(hobj,'ylim'))/dataratio(2));
    change = actualratio/availableratio;
    if change>1
        % we want a larger x/y ratio than given by the full axes
        % -> shrink y dimension
        pos0(2) = pos0(2) + psiz(2)*(1-1/change)/2;
        psiz(2) = psiz(2)/change;
    else
        % the contrary
        pos0(1) = pos0(1) + psiz(1)*(1-change)/2;
        psiz(1) = psiz(1)*change;
    end
    pos = [pos0 psiz];
end

 