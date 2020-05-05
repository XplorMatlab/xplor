classdef HeaderEdit < hgsetget
% function E = HeaderEdit(data [, callback])
% callback is a function with prototype fun(header) to execute once the ok
% button is pressed

properties
    % data (will be read-only)
    dat
    sz
    nd
    
    % header under construction
    cur_head = struct('sub_labels', cell(1,0), 'label', [], ...
        'unit', [], 'start', [], 'scale', [], 'values', [], ...
         'colors', [], 'isvalid', [], ...
        'confirmed', [], 'guessaction', '', 'all_guess', []);
    header

    % graphics
    hf
    table
    ok
    uconfirm
    contextmenu

    % callback
    callback
end

% Constructor
methods
function E = HeaderEdit(data, callback)
    % data
    if nargin < 1, data = rand(50, 40, 100); end
    if isa(data, 'xplr.Xdata')
        input_head = data.header;
        E.dat = data.data;
        E.sz = data.sz;
    else
        input_head = [];
        E.dat = data;
        E.sz = strictsize(data);
    end
    E.nd = length(E.sz);

    % callback
    if nargin >= 2
        E.callback = callback;
    end

    % guess header and store in a structure that admits invalid
    % definitions (units defined without start and scale being
    % defined)
    for i = 1:E.nd
        % guesses or other start value
        if isempty(input_head)
            candidates = xplr.Bank.get_recent_headers(E.sz(i), 20);
            okguess = ~isempty(candidates);
            if ~okguess
                candidates = xplr.Header('', E.sz(i));
            end
        else
            okguess = false;
            candidates = input_head(i);
        end
        % store in a different format more adequate for display
        clear candidates2
        for j = 1:length(candidates)
            head_j = candidates(j);
            sub_labels = {head_j.sub_labels.label};
            values = head_j.values;
            colors = [];
            kcolor = strcmp(sub_labels, 'ViewColor');
            if any(kcolor)
                colors = cell2mat(values(:, kcolor));
                sub_labels(kcolor) = [];
                values(:, kcolor) = [];
            end
            if isscalar(sub_labels) && strcmp(sub_labels{1}, head_j.label)
                % show sub_labels only if they are different from the
                % summary label
                sub_labels = [];
            end
            candidates2(j) = struct( ...
                    'sub_labels', {sub_labels}, 'label', head_j.label, ...
                    'unit', head_j.unit, 'start', head_j.start, 'scale', ...
                    head_j.scale, 'values', {values}, 'colors', colors, ...
                    'isvalid', true, 'confirmed', [], 'guessaction', ...
                    'reset', 'all_guess', [] ...
                ); %#ok<AGROW>
        end
        % current head value
        occurence = sum(E.sz(1:i) == E.sz(i));
        if length(candidates2) >= occurence, j = occurence; else j = 1; end
        E.cur_head(i) = candidates2(j);
        % setting data for guesses
        if okguess
            E.cur_head(i).confirmed = false;
            E.cur_head(i).guessaction = fn_switch( ...
                    isscalar(candidates2), 'reset', 'choose' ...
                );
            E.cur_head(i).all_guess = candidates2;
            % store all other alternatives as well
        else
            E.cur_head(i).confirmed = true;
        end
    end

    E.hf = figure('integerhandle', 'off', 'handlevisibility', 'off', ...
        'numbertitle', 'off', 'name', 'Set headers', ...
        'menubar', 'none', 'resize', 'off');
    set(E.hf, 'WindowButtonMotionFcn', @do_nothing)
    % force update of CurrentPoint when moving the mouse

    % menu for editing units
    menu_units(E)

    % init table
    init_table(E)
end
end

% Edit measures
methods
function menu_units(E)
    delete(findall( ...
            E.hf, 'type', 'uimenu', 'label', 'Edit recognized units' ...
        ))
    m = uimenu('parent', E.hf, 'label', 'Edit recognized units');
    measures = xplr.Bank.get_measures();
    labels = {measures.label};
    for i=1:length(labels)
        label = measures(i).label;
        units = {measures(i).units.unit};
        uimenu(m, 'label', [label ' (' fn_strcat(units, ',') ')'], ...
            'callback', @(u,e)edit_measure(E, label))
    end
    uimenu(m, 'label', 'Create a new set', 'separator', ...
        onoff(~isempty(labels)), 'callback', @(u,e)edit_measure(E, ''))
end

% this function is actually completely independent from E, except
% it re-displays the menus at the end
function edit_measure(E, old_label)

    % which measure to edit / new unit
    create_new = isempty(old_label);
    if create_new
        units = struct('unit', cell(1, 0), 'value', cell(1, 0));
    else
        measures = xplr.Bank.get_measures();
        idx = strcmp(old_label, {measures.label});
        units = measures(idx).units;
    end

    % figure
    hfm = figure('windowstyle', 'modal');
    clf(hfm, 'reset')
    set(hfm, 'numbertitle', 'off', 'name', 'Edit Units', ...
        'menubar', 'none', 'resize', 'off')

    % some constants for control positioning
    dx = 10;
    dy = 10;
    h = 14;
    htab = 250;
    W = 220;
    H = dy + h + 2*dy + htab + 2*dy + 2*h + dy;
    ycur = H;

    % use a fixed size for easy positioning of elements
    fn_setfigsize(hfm, [W, H]);
    set(hfm, 'defaultuicontrolunit', 'pixel')

    % label
    ycur = ycur - dy - h;
    uicontrol('style', 'text', 'string', 'Measure label', ...
        'position', [dx, ycur, (W - 2*dx)*.4, h]);
    x_label = uicontrol('style', 'edit', 'string', old_label, ...
        'position', [dx + (W - 2*dx)*.4, ycur, (W - 2*dx)*.6, h]);
    if create_new, set(x_label, 'ToolTipString', 'e.g. time'), end

    % data (add an empty line)
    data = [column({units.unit}), column({units.value})];

    % table
    ycur = ycur - 2*dy - htab;
    t = uitable('position', [dx, ycur, W - 2*dx, htab]);
    set( ...
        t, 'ColumnName', {'Unit', 'Value'}, 'ColumnFormat', ...
        {'char', 'numeric'}, 'ColumnEditable', true, 'RowName', '', ...
        'ColumnWidth', repmat({(W - 2*dx)/2 - 1}, 1, 2), 'Data', data ...
        )
    if create_new, set(t, 'ToolTipString', 'e.g. unit: ms, value: 1e-3'), end
    set(t, 'CellEditCallback', @checknewline), checknewline
    function checknewline(~, ~)
        data = get(t, 'Data');
        if isempty(data) || fn_find(data(end, :), 'any')
            data = [data; repmat({'', []}, 2, 1)];
            set(t, 'Data', data)
        end
    end

    % ok button
    ycur = ycur - 2*dy - 2*h;
    okmeas = uicontrol('string', 'ok', 'callback', @(u,e)delete(u), ...
        'position', [W - dx - 50, ycur, 50, 2*h]);

    % make units of controls 'normalized' for meaningful figure
    % resize behavior
    % set(get(hfm, 'children'), 'unit', 'normalized')

    % wait for use finished and grab data before closing window
    waitfor(okmeas)
    if ~ishandle(hfm), return, end % figure was closed, cancel
    label = get(x_label, 'String');
    data = get(t, 'Data');
    close(hfm)

    % save edited measure in bank
    while ~isempty(data) && ~fn_find(data(end, :), 'any')
        data(end, :)=[];
    end
    % remove empty lines

    if isempty(label) || isempty(data), return, end
    units = struct('unit', data(:, 1), 'value', data(:, 2));
    if create_new
        xplr.Bank.add_measure(label, units);
    else
        xplr.Bank.edit_measure(old_label, label, units);
    end

    % reinit menu to take into account the change in measures
    menu_units(E)
end
end

% Edit headers
methods
function init_table(E)
    % figure
    W = 560;
    H = 380;
    fn_setfigsize(E.hf, W, H)

    % ok button
    E.ok = uicontrol('parent', E.hf, 'string', 'ok', 'position', ...
        [W-80, 1, 80, 30], 'callback', @(u,e)done(E));

    % "confirm all" button masks ok button, but trigger done(E)
    % E.uconfirm = uicontrol('parent',E.hf,'string', ...
    % 'Confirm all','callback',@(u,e)confirm_all(E), ...
    % 'position',[W-80 1 80 30]);

    % reset button
    uicontrol('parent', E.hf, 'string', 'Reset all', 'callback', ...
        @(u,e)reset_all(E), 'position', [W-2*80, 1, 80, 30]);

    % empty table
    E.table = uitable('parent', E.hf, 'position', [1, 31, W, H-30], ...
        'ColumnName', { ...
            'Dim', 'Size', 'Label', 'Unit', 'Scale/Values', 'Colors', '' ...
            }, ...
        'ColumnFormat', { ...
            'numeric', 'numeric', 'char', 'char', 'char', 'char', 'char', ...
            }, ...
        'ColumnEditable', logical([0, 0, 1, 1, 1, 1, 0]));
    set(E.table, 'TooltipString', 'test')
    u = E.table;
    p = get(u, 'position');
    w = p(3);
    widths = {30, 55, 70, 70, [], 55, 55};
    w_avail = w - sum([widths{:}]) - 2;
    idx_auto = find(fn_isemptyc(widths));
    [widths{idx_auto}] = deal(floor(w_avail/length(idx_auto)));
    set(u, 'ColumnWidth', widths)
    set(u, 'RowName', [])
    set(u, 'CellEditCallback', @(u,e)celledit(E, e), ...
        'CellSelectionCallback', @(u,e)cell_select(E, e))

    % fill table
    display_header(E, 1:E.nd)
end
function display_header(E, idx)
    if nargin<2, idx = 1:E.nd; end
    t_data = get(E.table, 'Data');
    cur_data = t_data;
    [iL, iU, iV, iA, iC] = column_indices;
    for i = idx
        head = E.cur_head(i);
        % dimension number
        t_data{i, 1} = i;
        % length
        t_data{i, 2} = E.sz(i);
        % header info
        [t_data{i, [iL, iU, iV, iC]}] = display_headerinfo(head);
        % is guess?
        switch head.guessaction
            case ''
                error 'programming: guess action should be ''choose'' or ''reset'''
            case 'choose'
                t_data{i, iA} = 'Choose...';
            case 'reset'
                t_data{i, iA} = 'Reset';
        end
    end
    if ~isequal(t_data, cur_data), set(E.table, 'Data', t_data), end
    % which values are guessed and need being confirmed
    color_table(E)
end
function color_table(E)
    % color rows with values that need being confirmed
    col = zeros(E.nd, 3);
    allok = true;
    for i=1:E.nd
        if ~E.cur_head(i).isvalid
            col(i, :) = [1, .4, .4];
            allok = false;
        elseif ~E.cur_head(i).confirmed
            col(i, :) = [1, 1, 0];
            allok = false;
        elseif mod(i, 2)
            col(i, :) = [1, 1, 1];
        else
            col(i, :) = [1, 1, 1]*.94;
        end
    end
    set(E.table, 'BackgroundColor', col)
end
function celledit(E, e)
    i = e.Indices(1);
    cnames = get(E.table, 'ColumnName');
    ename = cnames{e.Indices(2)};

    % this automatically confirms the guess if there was a guess
    E.cur_head(i).confirmed = true;

    % update header
    switch ename
        case 'Label'
            update_label(E, i)
        case 'Unit'
            update_unit(E, i)
        case 'Scale/Values'
            update_value(E, i)
        case 'Colors'
            update_colors(E, i)
    end

    % update display
    check_valid(E, i)
    display_header(E, i)
end
function update_label(E, i)
    head = E.cur_head(i);
    t_data = get(E.table, 'Data');
    iL = column_indices;
    str = t_data{i, iL};
    tokens = regexp(str, '^(.*[^ ]) *\((.*)\)$', 'tokens');
    if ~isempty(tokens)
        % form 'label (label1*label2)'
        [head.label, sub_labels] = deal(tokens{1}{:});
        head.sub_labels = fn_strcut(sub_labels, '*');
    elseif any(str == '*')
        head.label = str;
        head.sub_labels = fn_strcut(str, '*');
    else
        head.label = str;
        head.sub_labels = [];
    end
    E.cur_head(i) = head;
end
function update_unit(E, i)
    t_data = get(E.table, 'Data');
    [~, iU] = column_indices;
    unit = t_data{i, iU};
    tokens = regexp(unit, '^(.*) \[(.*)\]$', 'tokens');
    if ~isempty(tokens), unit = tokens{1}{1}; end
    head = E.cur_head(i);
    % update start, scale, values if necessary
    if isempty(head.unit) && ~isempty(unit)
        [head.start, head.scale, head.values] = deal(1, 1, {});
    elseif ~isempty(head.unit) && isempty(unit)
        [head.start, head.scale, head.values] = deal([], [], {});
    end
    head.unit = unit;
    E.cur_head(i) = head;
end
function update_value(E, i)
    head = E.cur_head(i);
    t_data = get(E.table, 'Data');
    [~, ~, iV] = column_indices;
    [type, value] = read_value(t_data{i, iV}, E.sz(i));
    switch type
        case 'invalid'
            % do not accept change
        case 'measure'
            [head.scale, head.start] = dealc(value);
            head.values = [];
        case 'enum'
            [head.scale, head.start, head.values] = deal([]);
        case 'categorical'
            [head.scale, head.start] = deal([]);
            head.values = value;
    end
    E.cur_head(i) = head;
end
function update_colors(E, i)
    head = E.cur_head(i);
    t_data = get(E.table, 'Data');
    [~, ~, ~, ~, iC] = column_indices;
    % colors not allowed for 'measure' header
    if ~isempty(head.unit), return, end
    % no colors?
    if isempty(t_data{i, iC}), E.cur_head(i).colors = []; end
    % try setting colors
    try %#ok<TRYNC>
        colors = evalin('base', t_data{i, iC});
        if isequal(size(colors), [E.sz(i), 3])
            E.cur_head(i).colors = colors;
        end
    end
end
function check_valid(E, i)
    head = E.cur_head(i);
    % number of labels
    [n_label, n_column] = deal( ...
        length(head.sub_labels), size(head.values, 2) ...
        );
    if n_label == 0
        okv = (n_column <= 1);
    else
        okv = (n_column == n_label);
    end
    % if unit is defined, so must be start and scale
    okv = okv && (isempty(head.unit) ...
        || (~isempty(head.start) && ~isempty(head.scale)));
    E.cur_head(i).isvalid = okv;
end
function cell_select(E, e)
    % Special actions can occur when selecting a cell inside a row
    % where there is some guess
    [iL, iU, iV, iA, iC] = column_indices;
    if size(e.Indices, 1) ~= 1, return, end
    i = e.Indices(1);
    head_i = E.cur_head(i);
    if e.Indices(2) == iA
        % Perform 'guess action'
        switch head_i.guessaction
            case ''
                error 'programming: guess action should be ''choose'' or ''reset'''
            case 'choose'
                % Select among the list of all guesses: build and show
                % a menu with all possibilities
                deleteValid(E.contextmenu)
                m = uicontextmenu('parent', E.hf);
                E.contextmenu = m;
                n_guess = length(head_i.all_guess);
                for j=1:n_guess
                    try
                        [label, unit scale_value color] = ...
                            display_headerinfo(head_i.all_guess(j));
                        lab = fn_strcat( ...
                            {label, unit, scale_value color}, '; ' ...
                            );
                        uimenu(m, 'label', lab, 'callback', ...
                            @(u,e)use_guess(E, i, j))
                    catch
                        % probably a header was not saved properly
                        % in the bank, don't try to display it
                    end
                end
                uimenu(m, 'label', 'Reset', 'callback', ...
                    @(u,e)use_guess(E, i, 'reset'))
                set(m, 'position', get(E.hf, 'CurrentPoint'), ...
                    'visible', 'on')
            case 'reset'
                E.use_guess(i, 'reset')
        end
    elseif any(e.Indices(2) == [iL, iU, iV, iC]) && ~head_i.confirmed
        % Confirm guess
        head_i.confirmed = true; % line is confirmed in any case
        head_i.guessaction = 'choose';
        % allow user to go back to the guess
        E.cur_head(i) = head_i;
        % update display
        E.display_header(i)
    end
end
function use_guess(E, i, j)
    % update ore reset header
    all_guess = E.cur_head(i).all_guess;
    if strcmp(j, 'reset')
        head_i = E.cur_head(i);
        [ ...
            head_i.label, head_i.unit, head_i.start, head_i.scale, ...
            head_i.values, head_i.colors ...
            ] ...
            = deal('', '', [], [], cell(E.sz(i), 0), []);
    else
        head_i = all_guess(j);
    end
    % set guess data
    head_i.confirmed = true;
    % guess selected by user is automatically confirmed
    if isempty(all_guess) || ...
        (isscalar(all_guess) && ~strcmp(j, 'reset'))
        head_i.guessaction = 'reset';
    else
        head_i.guessaction = 'choose';
    end
    head_i.all_guess = all_guess;
    E.cur_head(i) = head_i;
    % update display
    display_header(E, i)
end
function confirm_all(E)
    [E.cur_head.confirmed] = deal(true);
    E.display_header()
end
function reset_all(E)
    for i = 1:E.nd
        E.use_guess(i, 'reset')
    end
end
function done(E)
    % build headers
    E.header = xplr.Header.empty(1,0);
    for i = 1:E.nd
        head = E.cur_head(i);
        if isempty(head.label)
            % set a label if there isn't
            head.label = ['dim', num2str(i)];
        end
        if isempty(head.scale)
            % categorical
            if isempty(head.sub_labels)
                head.sub_labels = {head.label};
            end
            if ~isempty(head.colors)
                head.sub_labels{end+1} = 'ViewColor';
                head.values(:, end+1) = num2cell(head.colors, 2);
            end
            if isempty(head.values)
                E.header(i) = xplr.Header(head.label, E.sz(i));
            else
                E.header(i) = xplr.Header(head.label, ...
                    head.sub_labels, head.values);
            end
        else
            % measure
            [unit, ~, measure, conversion] = read_unit(head.unit);
            if isempty(measure)
                dimlabel = xplr.Dimension_label( ...
                    head.label, 'numeric', unit ...
                    );
                conversion = 1;
            else
                dimlabel = xplr.Dimension_label( ...
                    head.label, 'numeric', measure.units ...
                    );
            end
            E.header(i) = xplr.Header( ...
                dimlabel, E.sz(i), head.start*conversion, ...
                head.scale*conversion ...
                );
        end
    end
    
    % close figure -> calling editHeader function can proceed
    delete(E.hf)
    drawnow
    
    % register headers to the bank
    if ~isempty(E.header)
        xplr.Bank.register_headers(E.header)
    end
    
    % execute callback if any
    if ~isempty(E.callback)
        E.callback(E.header);
    end
end
end
    
end

%---
function do_nothing(~, ~)
    % this function is set to the WindowButtonMotionFcn property of the figure
    % to force update of CurrentPoint when moving the cursor
end

%---
function [iL, iU, iV, iA, iC] = column_indices
    [iL, iU, iV, iA, iC] = deal(3,4,5,7,6);
end

%---
function [unit, str, measure, conversion] = read_unit(unit)
% if unit is recognized as a unit for a known measure, return an enhanced
% string display and the list of all units for this measure

    tokens = regexp(unit, '^(.*) \[(.*)\]$', 'tokens');
    if ~isempty(tokens), unit = tokens{1}{1}; end

    [measurelabel, conversion, measure] = xplr.Bank.get_unit_info(unit);
    isknown = ~isempty(measurelabel);
    if isknown
        comment = [' [' measurelabel ']'];
    else
        comment = [];
    end
    str = [unit, comment];

end

%---
function [type, value] = read_value(scale_value,n)

    % empty?
    if isempty(scale_value)
        type = 'enum';
        value = [];
        return
    end

    % scale + start?
    tokens = regexp(scale_value, '^(.*)\[ *(start){0,1}(.*)\]$', 'tokens');
    if ~isempty(tokens)
        tok = tokens{1};
        try %#ok<TRYNC>
            scale = evalin('base', tok{1});
            start = evalin('base', tok{3});
            type = 'measure';
            value = [scale, start];
            return
        end
    end

    % try evaluating in base workspace as a string that evaluates to a number
    try
        x = evalin('base', scale_value);
    catch
        x = [];
    end

    % number?
    if isscalar(x) && isnumeric(x) && n~=1
        type = 'measure';
        value = [x, 0];
        return
    elseif isvector(x) ...
        && isnumeric(x) ...
        && length(x) == n ...
        && max(abs(diff(x, 2))) < abs(diff(x(1:2)))/1e6
        % vector of equally spaced values
        if max(abs(diff(x, 2))) == 0
            scale = diff(x(1:2));
        else
            scale = (x(n)-x(1))/(n-1);
        end
        start = x(1);
        type = 'measure';
        value = [scale, start];
        return
    end

    % list of items?
    if isvector(x) && n > 1, x = column(x); end
    if isnumeric(x) && size(x, 1) == n
        list = num2cell(x);
    elseif iscell(x) && size(x, 1) == n
        list = x;
    else
        sep = fn_switch(any(scale_value == ','), ',', ' ');
        list = fn_strcut(scale_value, sep);
        list = fn_regexptokens(list, '^ *(.*?) *$');
        % remove blanks at the beginning and end
        if isvector(list) && n>1, list = column(list); end
    end
    if length(list)==n
        type = 'categorical';
        value = list;
        return
    else
        % failed interpreting scale_value
        type = 'invalid';
        value = [];
        return
    end

end

%---
function [label, unit, scale_value, color] = display_headerinfo(head)

    % label
    if isempty(head.sub_labels)
        label = head.label;
    else
        autolabel = fn_strcat(head.sub_labels, '*');
        if strcmp(autolabel, head.label)
            label = head.label;
        else
            label = [head.label, ' (', autolabel, ')'];
        end
    end
    % unit
    [~, unit] = read_unit(head.unit);
    % scale/values
    if isempty(head.scale)
        scale_value = display_value('categorical', head.values);
    else
        scale_value = display_value('measure', [head.scale, head.start]);
    end
    % color
    if isempty(head.colors)
        color = '';
    else
        color = sprintf('[%ix%i array]', size(head.colors));
    end

end

%---
function str = display_value(type, value)

switch type
    case 'measure'
        str = [num2str(value(1), 12), ' [start ', num2str(value(2),12), ']'];
    case 'categorical'
        if isempty(value)
            str = '';
        elseif isvector(value)
            str = fn_strcat(value, ',');
        else
            str = sprintf('%ix%i table', size(value));
        end
end

end
