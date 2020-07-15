function el = pixelposlistener(source,varargin)
% function el = pixelposlistener(source,[target,],callback)
%---
% Add a listener that will execute whenever the pixel position of an object
% is changed. 
% In Matlab version R2014b and later, this just adds a listener to the
% object 'LocationChanged' event. In earlier versions, this is a wrapper
% for pixelposwatcher class.
%
% See also brick.pixelpos, brick.pixelsizelistener

% Thomas Deneux
% Copyright 2015-2017

% Input
switch nargin
    case 2
        callback = varargin{1};
        target = [];
    case 3
        [target callback] = deal(varargin{:});
end

% Create listener
el = [brick.connect_listener(source,target,'LocationChanged',callback) brick.connect_listener(source,target,'SizeChanged',callback)];

% Output?
if nargout==0, clear el, end
