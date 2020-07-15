function el = pixelsizelistener(source,varargin)
% function el = pixelsizelistener(source,[target,],callback)
%---
% Add a listener that will execute whenever the pixel size of an object
% is changed.
%
% See also brick.pixelsize, brick.pixelposlistener

% Thomas Deneux
% Copyright 2015-2020

% Input
switch nargin
    case 2
        callback = varargin{1};
        target = [];
    case 3
        [target, callback] = deal(varargin{:});
end

% Create listener
el = brick.connect_listener(source,target,'SizeChanged',callback);

% Output?
if nargout==0, clear el, end
