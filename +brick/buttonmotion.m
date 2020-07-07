function varargout = buttonmotion(fun,varargin)
%BUTTONMOTION Execute a task while mouse pointer is moved (try 'brick.buttonmotion demo')
%---
% function varargout = buttonmotion(fun[,hf][,'doup'][,'pointer',pointer])
% function [moved ...] = buttonmotion(...,'moved?')
% function buttonmotion('demo')
%---
% utility for executing a task while mouse button is pressed
% 
% See also brick.moveobject

% Thomas Deneux
% Copyright 2007-2017

% Note on the implementation: unfortunately, it is not possible to simply
% use the figure 'BusyAction' and 'Interruptible' properties. Indeed, they
% act on all callbacks, whereas we would like to treat differently the
% motion callbacks (cancel if system is busy) and the button press
% callbacks (queue if system is busy).

if nargin==0, help brick.buttonmotion, return, end


disp_if_debug('entering brick.buttonmotion function')

% Input
if ischar(fun) && strcmp(fun,'demo')
    demo()
    return
end
hf = []; doup = false; pointer = ''; checkmove = false;
k = 0;
while k<length(varargin)
    k = k+1;
    a = varargin{k};
    if ischar(a)
        switch a
            case 'doup'
                doup = true;
            case 'pointer'
                k = k+1;
                pointer = varargin{k};
            case 'moved?'
                checkmove = true;
            otherwise
                error 'unknown flag'
        end
    else
        hf = a;
    end
end
if isempty(hf), hf = gcf; end

% How many outputs of the function are expected
nout = max(0,nargout-checkmove);

% brick.buttonmotion did not terminate correctly last time?
motionfcn = get(hf,'WindowButtonMotionFcn');
if iscell(motionfcn) && isequal(motionfcn{1},@motionexec)
    % check whether there was a queued stop command
    stopdebugstr = getappdata(hf,'brick_buttonmotion_queuestop');
    % terminate previous call in any case
    disp_if_debug('New call to brick.buttonmotion in this figure before the previous one was finished: terminate the previous one first')
    terminate(hf)
    disp_if_debug('Now go on with the new call to brick.buttonmotion')
else
    stopdebugstr = '';
end

% Backup properties that will be overriden
state = brick.get(hf,{'WindowButtonDownFcn' 'WindowButtonUpFcn' 'WindowButtonMotionFcn' 'Interruptible' 'Pointer'});
setappdata(hf,'brick_buttonmotion_savestate',state)

% Motion
set(hf,'WindowButtonMotionFcn',{@motionexec 'move' fun nout}, ...
    'WindowButtonUpFcn',{@motionexec 'stop'},'WindowButtonDownFcn',{@motionexec 'stop2'})
set(hf,'Interruptible','on') % make sure events will be interruptible, to not miss the 'stop' event
if ~isempty(pointer), set(hf,'pointer',pointer), end
setappdata(hf,'brick_buttonmotion_scrolling',true)
setappdata(hf,'brick_buttonmotion_busy',false)
setappdata(hf,'brick_buttonmotion_moved',false)
setappdata(hf,'brick_buttonmotion_lastmoverejected',false)
if nout, setappdata(hf,'brick_buttonmotion_output',cell(1,nout)), end

% Wait for motion end
if ~isempty(stopdebugstr)
    % if there was a queued stop command do not start
    disp_if_debug(['The previous call had queued command ' stopdebugstr ', so stop new call immediately'])
else
    disp_if_debug('waiting for motion end')
    waitfor(hf,'WindowButtonMotionFcn') % at motion end, function terminate will be executed, and this property will be changed
    disp_if_debug('finished waiting')
end
    
% Finish
if ~ishandle(hf), return, end % figure has been closed in the mean while
if doup || any(getappdata(hf,'brick_buttonmotion_lastmoverejected'))
    setappdata(hf,'brick_buttonmotion_scrolling',false)
    exec(fun,nout,hf);
    try rmappdata(hf,'brick_buttonmotion_scrolling'), end %#ok<TRYNC>
end
try rmappdata(hf,'brick_buttonmotion_lastmoverejected'), end %#ok<TRYNC>
if checkmove
    varargout = {getappdata(hf,'brick_buttonmotion_moved')};
else 
    varargout = {}; 
end
try rmappdata(hf,'brick_buttonmotion_moved'), end %#ok<TRYNC>
if nout
    varargout = [varargout getappdata(hf,'brick_buttonmotion_output')];  %#ok<VARARG>
    try rmappdata(hf,'brick_buttonmotion_output'), end %#ok<TRYNC>
end 

%---
function motionexec(hf,~,actionflag,fun,nout)

persistent kid
if isempty(kid), kid = 0; end
kid = brick.mod(kid+1,1000);

% start execution
debugstr = [actionflag ' ' num2str(kid)];
disp_if_debug(['start ' debugstr])
if strcmp(actionflag,'stop2'), actionflag = 'stop'; end

% custom queuing/canceling system
if getappdata(hf,'brick_buttonmotion_busy')
    switch actionflag
        case 'move' % cancel
            setappdata(hf,'brick_buttonmotion_lastmoverejected',true)
            disp_if_debug(['rejct ' debugstr])
            return
        case 'stop' % queue
            setappdata(hf,'brick_buttonmotion_queuestop',debugstr)
            disp_if_debug(['queue ' debugstr])
            return
    end
elseif strcmp(actionflag,'move')
    setappdata(hf,'brick_buttonmotion_lastmoverejected',false)
end
setappdata(hf,'brick_buttonmotion_busy',true)


% stop (not queued)
if strcmp(actionflag,'stop')
    disp_if_debug(['end   ' debugstr])
    terminate(hf)
    return
end


% evaluate function
disp_if_debug(['exec  ' debugstr])
try 
    exec(fun,nout,hf); 
catch ME
    terminate(hf)
    rethrow(ME)
end


% end of current execution
disp_if_debug(['end   ' debugstr])
drawnow % allow queued events to be processed (and canceled, because of 'brick.buttonmotion_busy' flag)
setappdata(hf,'brick_buttonmotion_moved',true)
setappdata(hf,'brick_buttonmotion_busy',false)

% stop (queued)
try
    debugstr = getappdata(hf,'brick_buttonmotion_queuestop');
    if ~isempty(debugstr)
        disp_if_debug(['exec  queued ' debugstr])
        terminate(hf)
    end
end

%---
function exec(fun,nout,hf)

if nargin>=2 && nout, out = cell(1,nout); end
if ischar(fun)
    evalin('base',fun);
elseif isa(fun,'function_handle')
    if nout, [out{:}] = feval(fun); else feval(fun); end
elseif iscell(fun)
    if nout, [out{:}] = feval(fun{:}); else feval(fun{:}); end
else
    error bad
end
if nout, setappdata(hf,'brick_buttonmotion_output',out), end


%---
function terminate(hf)

state = getappdata(hf,'brick_buttonmotion_savestate');
brick.set(hf,state)
% remove data attached to the figure, except 'moved', 'lastmoverejected'
% and 'output' which will be needed in the main brick.buttonmotion after
% terminate(hf) has executed
rmappdata(hf,'brick_buttonmotion_savestate')
rmappdata(hf,'brick_buttonmotion_busy')
if isappdata(hf,'brick.buttonmotion_queuestop'), rmappdata(hf,'brick_buttonmotion_queuestop'), end
rmappdata(hf,'brick_buttonmotion_scrolling')

%---
function disp_if_debug(varargin)

if eval('false')
    str = [];
    for k=1:nargin
        x = varargin{k};
        if iscell(x)
            if isa(x{1},'function_handle')
                x = func2str(x{1});
            else
                error('don''t know how to display cell array')
            end
        end
        str = [str x]; %#ok<AGROW>
    end
    disp(str)
end

%---
function demo

C = {'figure(1), clf'
    'ht = uicontrol(''style'',''text'',''backgroundcolor'',''y'');'
    'fun = @()set(ht,''string'',num2str(get(1,''CurrentPoint'')));'
    'set(1,''buttondownfcn'',@(u,e)brick.buttonmotion(fun,''doup''))'};
brick.dispandexec(C)


