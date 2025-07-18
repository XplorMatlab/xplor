function out = progress(varargin)
% function progress(prompt,max[,ht][,'ignoresub'][,'noerase'])
% function progress(prompt,'%',[,ht][,'ignoresub'][,'noerase'])
% function progress(text)
% function progress(i)
% function progress(i,'pause')
% function progress('end')
% function progress('cont')
% function progress('in',ht|[])
% function progress('screen')
% function progress('elapsed|elapsedmin')
% function x = progress(prompt,...)
% function progress(i,x)
%---
% progress indicator
%
% Examples:
%   - Simple call to brick.progress
%
%     n_item = 5
%     brick.progress('item', n_item)
%     for i=1:n_item
%         brick.progress(i)
%         pause(.25) % do some stuff
%     end
%
%   - Nested calls to brick.progress
%
%     n_items = [3, 4, 2];
%     n_batch = length(n_items);
%     brick.progress('batch', n_batch)
%     for i = 1:n_batch
%         brick.progress(i)
%         n_item_i = n_items(i);
%         brick.progress('item', n_item_i)
%         for j=1:n_item_i
%             brick.progress(j)
%             pause(.25) % do some stuff
%         end
%         brick.progress('end') % needed to signal that we finished
%     end
%
% See also pg

% Thomas Deneux
% Copyright 2003-2017 

persistent x            % structure with persistent information
persistent ht0          % default place for displaying progress (handle or '' for command prompt)

% First check fast whether too few time ellapsed since last display and we
% should not loose time making another display
if ~ischar(varargin{1})
    i = varargin{1};
    do_pause = (nargin>1 && strcmp(varargin{2},'pause'));
    if toc(x(1).lastdisp)<.1 && i~=x(1).max && ~do_pause
        return
    end
end

% Detect nested call to brick.progress (variable calllevel states the level
% of nesting)
stack = dbstack; 
if isscalar(stack), caller = ''; else caller = stack(2).name; end
if nargout==1
    % syntax x = brick.progress(prompt,...)
    calllevel = 0;
elseif nargin>=2 && isstruct(varargin{2})
    % syntax brick.progress(i,x)
    calllevel = 1;
    x = varargin{2};
elseif isempty(x) || ~isfield(x,'caller')
    calllevel = 0; 
elseif strcmp(caller,x(1).caller) 
    calllevel = 1; 
elseif length(x)>1 && strcmp(caller,x(2).caller)
    calllevel = 2; %#ok<NASGU> [call level is 2, by doing x(1)=[] we will make it 1]
    % we got back one level up, which means we finished with the most nested call!
    if x(1).ignoresub
        disp 'it seems that a initializing call to brick.progress with flag ''ignoresub'' was not properly terminated with a call ''brick.progress(''end'')'''
    end
    x(1) = [];
    x(1).doerase = false; % since there was some printing from a sub-level, we don't want the super-level to erase anymore
    calllevel = 1;
else
    calllevel = 0; % this is either a new call, or we lost track of the right parameters
    if x(1).ignoresub
        % make sure that the registered caller is in the stack
        if any(strcmp({stack.name},x(1).caller))
            % good
            return
        else
            disp 'entry to brick.progress with flag ''ignoresub'' was not completed with ''brick.progress('end')'''
            x(1) = [];
            if length(x)>1 && strcmp(caller,x(2).caller)
                x(1).doerase = false; % since there was some printing from a sub-level, we don't want the super-level to erase anymore
                calllevel = 1;
            end
        end
    end
end


if ischar(varargin{1})
    prompt = varargin{1};
    
    % SPECIAL COMMANDS
    if nargin==1
        switch prompt
            case 'end'
                if isempty(x), return, end
                if ishandle(x(1).ht)
                    set(x(1).ht,'string','')
                else
                    fprintf(repmat('\b',1,x(1).promptsize+1+x(1).isize+length(x(1).after)+1))
                end
                x(1) = [];
                return
            case 'cont'
                fprintf(repmat(' ',1,x(1).promptsize+1+x(1).isize+length(x(1).after)))
                return
            case 'elapsed'
                disp(['elapsed ' num2str(toc(x(1).timer)) 's.'])
                return
            case 'elapsedmin'
                disp(['elapsed ' num2str(toc(x(1).timer)/60,'%.1f') 's.'])
                return
            case 'screen'
                if ~isempty(x) && any(ishandle(x(1).ht)), set(x(1).ht,'string',''), end
                ht0 = [];
                return
        end
    elseif strcmp(prompt,'in')
        ht0 = varargin{2};
        return
    end
    
    % SINGLE TEXT DISPLAY
    if nargin<2
        if ~isempty(x), [x.doerase] = deal(false); end % since we are printing something, we don't want to erase anymore
        updatedisplay(-1,varargin{1},ht0)
        return
    end
    
    % INITIALIZATION
    % set current structure and remember last parameters in case of nested calls to brick.progress
    x0 = struct( ...
        'prompt',[], ...       % prompt 
        'promptsize',[], ...   % size of prompt - used for 'end' and 'cont'
        'curi',[], ...         % remembers current progress - avoids printing again the same thing
        'after',[], ...        % message after the progress - e.g. %, or /139
        'isize',[], ...        % length of the string for progress number
        'format',[], ...       % format for the number
        'doerase',[], ...      % do erase previous progress - on by default, at init use negative 'max' number to set it off)
        'pflag',[], ...        % use % rather than /139
        'max',[], ...          % maximal progress number
        'ht',[], ...           % current place for displaying progress - if specified at init, is used instead of ht0
        'timer',[], ...        % timer started at init - used for 'elapsed' and 'elapsedmin'
        'lastdisp',[], ...     % timer started at last display
        'caller',[], ...       % name of caller function ('' for base workspace) - used to prevent nested calls to progress
        'ignoresub',[] ...     % ignore subsequent calls to brick.progress from another caller until brick.progress('end') has been invoked
        );
    if isempty(x), x = x0; else x = [x0 x(1)]; end
    % use caller to detect nested calls
    stack = dbstack;
    if isscalar(stack), x(1).caller = ''; else x(1).caller = stack(2).name; end
    % sizes
    x(1).promptsize = length(prompt);
    i = 0; x(1).curi = i;
    x(1).pflag = strcmp(varargin{2},'%');
    if x(1).pflag
        x(1).max = 1;
        x(1).isize = 3;
    else
        x(1).max = varargin{2};
        if ischar(x(1).max)
            x(1).max = str2double(x(1).max);
        end
        x(1).max = double(abs(x(1).max));
        x(1).isize = floor(log10(x(1).max))+1;
    end
    if ~isempty(ht0) && ~ishandle(ht0), ht0=[]; end
    x(1).ht = ht0;
    x(1).ignoresub = false;
    for karg = 3:length(varargin)
        a = varargin{karg};
        if ischar(a) 
            switch a
                case 'ignoresub'
                    x(1).ignoresub = true;
                case 'noerase'
                    x(1).doerase = false;
            end
        else
            x(1).ht = a;
        end
    end
    if isempty(x(1).doerase), x(1).doerase = isempty(x(1).ht) && (x(1).max>0); end
    % format
    x(1).prompt = prompt;
    x(1).format = ['%' num2str(x(1).isize) 'i'];
    if x(1).pflag
        x(1).after = '%%';
    else
        x(1).after = ['/' sprintf(x(1).format,x(1).max)];
    end
    % initial display
    updatedisplay(-1,[prompt ' ' sprintf(x(1).format,0) x(1).after],x(1).ht)
    % start timers
    x(1).timer = tic;
    x(1).lastdisp = tic-1e6; % trick: save the time it was 1s before, so that next display attempt will detect a delay > 1s and will not cancel the display
    % output?
    if nargout==1
        out = x;
    end
else   
    % STATE PROGRESS
    % input
    i = varargin{1};
    do_pause = (nargin>1 && strcmp(varargin{2},'pause'));
    % did we lose track of the right parameters?
    if ~calllevel
        if isscalar(x), disp 'strange: unregistered caller', end 
        str = ['???? ' num2str(i) '/??'];
        updatedisplay(-1,str,'')
        return
    end
    % display
    if x(1).pflag, i = floor(i/x(1).max*100); end
    if (i == x(1).curi), return, end
    x(1).curi = i;
    if x(1).doerase
        nerase = x(1).isize+brick.switch_case(x(1).pflag,1,length(x(1).after));
        updatedisplay(nerase,[sprintf(x(1).format,i) x(1).after],x(1).ht)
    else
        updatedisplay(-1,[x(1).prompt ' ' sprintf(x(1).format,i) x(1).after],x(1).ht)
    end
    % pause after display
    if do_pause, pause, end
    % reinit timer
    x(1).lastdisp = tic;
end

drawnow

%---
function updatedisplay(nerase,str,ht)

if ~isempty(ht) && ishandle(ht)
    if nerase>0
        s = get(ht,'string');
        ns = length(s);
        s = s(1:ns-nerase);
    else
        s = [];
    end
    set(ht,'string',[s str])
    drawnow
else
    fprintf([repmat('\b',1,nerase+1) str '\n'])
end
