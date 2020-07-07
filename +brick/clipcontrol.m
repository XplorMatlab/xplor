function hu = clipcontrol(ha)
%CLIPCONTROL A wrapper of brick.sensor that controls the clipping range applied to an image
%---
% function hu = clipcontrol(ha)
%---
% This is a wrapper for brick.sensor that controls clipping of axes ha (ha can
% be a vector of several axes handle)

% Thomas Deneux
% Copyright 2015-2017

if nargin<1, ha = gca; end

if ~isempty(get(ha(1),'deletefcn')), error 'axes already has a delete function', end


hf = get(ha(1),'parent');
clip = get(ha(1),'clim');
set(ha,'clim',clip)
hu = brick.sensor('value',clip,'callback',@(u,e)set(ha,'clim',u.value));
brick.controlpositions(hu,hf,[0 1 0 0],[5 -20 100 15])
set(ha(1),'deletefcn',@(u,e)delete(hu(ishandle(hu))))
