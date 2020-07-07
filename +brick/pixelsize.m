function [out1 out2] = pixelsize(hobj,varargin)
% function siz = pixelsize(hobj[,'strict'])
% function [w h] = pixelsize(hobj[,'strict'])
%---
% returns the width and height in pixels of any object without needing to
% change any units values
%
% In R2014b and later, wraps function getpixelposition
%
% See also brick.pixelpos, brick.pixelposlistener, brick.pixelsizelistener,
% brick.objectsize

% Thomas Deneux
% Copyright 2011-2019

% strict?
strict = brick.flags({'strict'},varargin);

if strict
    pos = brick.pixelpos(hobj,'strict');
    siz = pos(3:4);
else
    pos = getpixelposition(hobj);
    siz = pos(3:4);
end

if nargout==2
    out1 = siz(1);
    out2 = siz(2);
else
    out1 = siz;
end