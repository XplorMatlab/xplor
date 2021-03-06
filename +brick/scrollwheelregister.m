function scrollwheelregister(hobj,arg,arg2)
%SCROLLWHEELREGISTER Define scrollwheel actions specific to which object the mouse is hovering over
%---
% function scrollwheelregister(hobj,callback[,'default|off'])
% function scrollwheelregister(hobj,'')
% function scrollwheelregister(hobj,'off|on')
% function scrollwheelregister(hobj,'default')
% function scrollwheelregister(hobj,'mask')
% function scrollwheelregister(hobj,'repair')
%---
% This function can register different scroll wheel actions for different
% objects children of the same figure, even though Matlab provides only one
% 'WindowScrollWheelFcn' callback function for the whole figure. Which
% action will be executed is determined on the position of the mouse in the
% figure.
%
% Input:
% - hobj        handle of figure, axes, uipanel or uicontrol
% - callback    function with prototype fun(n), where n is the number of
%               vertical scrolls; use empty string to unregister callback,
%               of fun(n,modifiers), where modifier is a cell array of the
%               modifier keys that are pressed during the scrolling
% - 'default'   make the scroll wheel action associated to this object the
%               figure default action
% - 'on|off'    inactivate or reactivate the callback associated to object
%               hobj
% - 'mask'      the object has not associated callback, but the default
%               scroll wheel action will not be executed when the mouse is
%               over this object
% - 'repair'    repaire registration mechanism that might have been damaged
%               consecutive to an error, or to setting of window
%               scrollwheelfcn property by another function
%
% See also windowcallbackmanager

% Thomas Deneux
% Copyright 2012-2017

if nargin==0, help brick.scrollwheelregister, return, end

% parent figure
hf = brick.parentfigure(hobj);

% window callback manager for this figure
W = brick.windowcallbackmanager(hf);

% execute the appropriate method
if ischar(arg)
    switch arg
        case ''
            unregister(W,hobj)
            return
        case {'on' 'off'}
            setactive(W,hobj,brick.boolean(arg))
            return
        case 'default'
            setdefault(W,hobj)
            return
        case 'mask'
            register(W,hobj,[])
            return
        case 'repair'
            repair(W)
            return
    end
end
callback = arg;
active = (nargin<3 || ~strcmp(arg2,'off'));
register(W,hobj,callback,active)
if nargin>=3 && strcmp(arg2,'default')
    setdefault(W,hobj)
end
