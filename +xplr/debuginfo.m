function out = debuginfo(category, message, varargin)
%DEBUGINFO  display debug messages for developers only
%---
% function debugactive = debuginfo()
% function debuginfo(code, message [, sprintf arguments])
% function debuginfo('select')
% function debuginfo('on|off')
% function debuginfo('reset')
%---
% This function is only for XPLOR developers.
% First activate debug info display by typing debuginfo('on') in Matlab.
% Each message belongs to a category, and categories to display or not can
% be selected by the user by typing debuginfo('select').
% Each time debuginfo(category, message) is called with a new category
% name, user is asked whether this category should be displayed.
% If some register categories are not used any more, use debuginfo('reset')
% to remove all categories and rebuild the list from upcoming call to
% debuginfo(category, message).
%
% Special category 'stop' stops the debugger in addition to displaying the
% message.

switch nargin
    case 0
        out = debug_active();
    case 1
        flag = category;
        switch(flag)
            case {'on' 'off'}
                set_active(strcmp(flag,'on'))
            case 'select'
                select_categories()
            otherwise
                error('unknown flag ''%s''', flag)
        end
    otherwise
        display_info(category, message, varargin{:})
end

%---
function state = get_state(set_value)
% keep state in memory to avoid multiple reading from file

persistent memstate

if nargin==1
    memstate = set_value;
    return
end

if isempty(memstate)
    memstate = fn_userconfig('xplr.debuginfo');
    if isempty(memstate)
        memstate = struct('active', false, 'categories', struct());
        save_state(memstate)
    end
end

state = memstate;


%---
function save_state(state)

% save on disk
fn_userconfig('xplr.debuginfo',state)

% replace value in get_state persistent variable
get_state(state)


%---
function ok = debug_active()

state = get_state();
ok = state.active;


%---
function display_info(category, message, varargin)

state = get_state();
if ~state.active, return, end

try
    do_display = state.categories.(category);
catch
    % first time this category is met
    quest = sprintf('Do you want to display debug information of new category ''%s''?', ...
        category);
    do_display = fn_dialog_questandmem(quest,'xplr.debuginfo');
    
    % memorize answer
    state.categories.(category) = do_display;
    save_state(state)
end

% display if appropriate
if do_display
    fprintf([category ': ' message '\n'], varargin{:})
    
    % special: stop the debugger if category is 'stop'
    if strcmp(category, 'stop')
        keyboard
    end
end


%---
function set_active(value)

state = get_state();
state.active = value;
save_state(state)
if value
    disp 'xplor debug information enabled'
    disp 'select which categories to display with xplr.debuginfo(''select'')'
else
    disp 'xplor debug information disabled'
end


%---
function select_categories()

state = get_state();
names = fieldnames(state.categories);
selected = struct2array(state.categories);
if isempty(names)
    disp('xplr.debuginfo: no message category met yet')
    return
end
    
state =get_state();
[selection, ok] = listdlg( ...
    'Name', 'xplr.debuginfo', ...
    'PromptString', 'Select which debuginfo categories should be displayed', ...
    'ListString', names, ...
    'InitialValue', find(selected));

if ok
    selected(:) = false;
    selected(selection) = true;
    state.categories = cell2struct(num2cell(selected), names, 2);
    save_state(state)
end

        
    


