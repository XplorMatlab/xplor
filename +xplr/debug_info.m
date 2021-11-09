function out = debug_info(category, message, varargin)
%debug_info  display debug messages for developers only
%---
% function debugactive = debug_info()
% function debug_info(code, message [, sprintf arguments])
% function debug_info('select')
% function debug_info('on|off|file')
% function debug_info('reset')
%---
% This function is only for XPLOR developers.
% First activate debug info display by typing debug_info('on') in Matlab.
% Each message belongs to a category, and categories to display or not can
% be selected by the user by typing debug_info('select').
% Each time debug_info(category, message) is called with a new category
% name, user is asked whether this category should be displayed.
% If some register categories are not used any more, use debug_info('reset')
% to remove all categories and rebuild the list from upcoming call to
% debug_info(category, message).
%
% Special category 'stop' stops the debugger in addition to displaying the
% message.

switch nargin
    case 0
        out = debug_active();
    case 1
        flag = category;
        switch(flag)
            case {'on', 'off', 'file'}
                set_active(flag)
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

persistent mem_state

if nargin==1
    mem_state = set_value;
    return
end

if isempty(mem_state)
    mem_state = brick.userconfig('xplr.debug_info');
    if isempty(mem_state)
        mem_state = struct('active', 'off', 'categories', struct());
        save_state(mem_state)
    end
end

state = mem_state;


%---
function save_state(state)

% save on disk
brick.userconfig('xplr.debug_info', state)

% replace value in get_state persistent variable
get_state(state)


%---
function flag = debug_active()

state = get_state();
flag = state.active;


%---
function display_info(category, message, varargin)

state = get_state();
switch state.active
    case 'off'
        return
        
    case 'on'
        % check whether category is displayed
        try
            do_display = state.categories.(category);
        catch
            % first time this category is met
            quest = sprintf('Do you want to display debug information of new category ''%s''?', ...
                category);
            do_display = brick.dialog_questandmem(quest,'xplr.debug_info');

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

    case 'file'
        xplr.log_to_file([category ': ' message])
        
end

%---
function set_active(value)

state = get_state();
state.active = value;
save_state(state)
switch value
    case 'on'
        disp 'xplor debug information enabled'
        disp 'select which categories to display with xplr.debug_info(''select'')'
    case 'off'
        disp 'xplor debug information disabled'
    case 'file'
        disp 'all xplor debug information will be logged to file xplor.log'
end


%---
function select_categories()

state = get_state();
names = fieldnames(state.categories);
selected = cell2mat(struct2cell(state.categories));
if isempty(names)
    disp('xplr.debug_info: no message category met yet')
    return
end
    
state = get_state();
[selection, ok] = listdlg( ...
    'Name', 'xplr.debug_info', ...
    'PromptString', 'Select which debug_info categories should be displayed', ...
    'ListString', names, 'SelectionMode', 'multiple', ...
    'CancelString', 'None', ...
    'InitialValue', find(selected));

selected(:) = false;
if ok, selected(selection) = true; end
state.categories = cell2struct(num2cell(selected), names);
save_state(state)
