function hu=okbutton(varargin)
%OKBUTTON Small 'ok' button waits to be pressed
%---
% function okbutton([properties...])
% function okbutton('wait'[, properties...])
% function okbutton('closefigure'[, properties...])
%---
% Create a small 'ok' button at the bottom-left corner of the current
% figure (use 'parent' property to specify a specific figure), which will
% auto-delete when being pressed.
%
% Example usage:
%     figure
%     imdata = load('clown.mat');
%     imshow(imdata.X,imdata.map)
%     title 'Draw rectangle sub-region to look at. Press ok when done.'
%     ok = brick.okbutton();
%     hl = drawrectangle();       % returns once a rectangle is drawn, but user can further modify the rectangle
%     waitfor(ok)                 % wait until ok button has been pressed
%     rect = get(hl, 'Position'); % get the coordinates of the rectangle
%     close(gcf) 
%
% See also brick.nextbutton

% Thomas Deneux
% Copyright 2009-2017

if nargin>=1
    flag = varargin{1};
    if brick.ismemberstr(flag,{'wait' 'closefigure'})
        varargin(1) = [];
    else
        flag = '';
    end
else
    flag = '';
end

hu = uicontrol('position',[10 10 30 15],'string','ok', ...
    'callback','delete(gcbo)',varargin{:});
switch flag
    case 'wait'
        waitfor(hu)
    case 'closefigure'
        set(hu,'callback','close(gcf)')
    case ''
    otherwise
        error('unknown flag ''%s''',flag)
end
if nargout==0, clear hu, end
