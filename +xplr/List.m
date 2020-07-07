classdef List < xplr.GraphNode
    % function L = list(filter[,'in',uipanel][,other options...])
        
    properties (SetAccess='private')
        F           % xplr.filterAndPoint object
        sel_type
        h_list
        h_p  % direct parent
        h_f  % figure parent
        h_label
        menu
        value_str    % precomputed list of values
    end
    properties (SetObservable, AbortSet=true)
        sel_mult_in = true;
        scroll_wheel = 'on'; % 'on', 'off' or 'default'
        selection_prompt_name = 'none'; % 'all', 'groups' or 'none'
    end
    properties
        selection_save_file       % name of file for saving current selection
    end

    % Constructor and Destructor
    methods
        function L = List(F, varargin)
            % options for initialization
            opt = struct( ...
                'in',                   [] ...
                );
            if nargin == 0
                head = xplr.Header('testheader', 10);
                F = xplr.FilterAndPoint(head);
            end
            [opt, opt_add] = parse_input(opt, varargin{:});
            
            % check filter
            L.F = F;
            if F.nd_in ~= 1 || F.nd_out ~= 1, error 'input and output of list filter must be one-dimensional', end
            if ~isa(F,'xplr.FilterAndPoint'), error 'list can act only on a filterAndPoint object', end
            L.sel_type = brick.switch_case(F.header_in.categorical, 'indices', 'point1D');
            
            % watch filter deletion
            addlistener(L.F, 'ObjectBeingDestroyed', @(u,e)delete(L));
            
            % add 'soft_selection' label to output header
            F.augment_header('soft_selection', 'logical')
            
            % uipanel container
            new_figure = isempty(opt.in);
            if new_figure
                % create uipanel in a new figure
                L.h_p = uipanel('parent', figure);
            elseif strcmp(get(opt.in, 'type'), 'uipanel')
                L.h_p = opt.in;
            else
                error 'input container must be an uipanel object'
            end
            
            % parent figure
            L.h_f = brick.parentfigure(L.h_p);

            % create several components
            % (list)
            L.h_list = uicontrol('parent', L.h_p, 'style', 'listbox', 'min', 0, 'max', 2, ...
                'callback', @(h_list, evnt)event(L, 'select'), ...
            	'keypressfcn', @(h_list, evnt)keypress(L, evnt));
            brick.controlpositions(L.h_list, L.h_p, [0, 0, 1, 1], [8, 5, -16, -5-21-2])
            % (label)
            L.h_label = uicontrol('parent', L.h_p, 'style', 'text', ...
                'string', L.F.header_out.label, ...
                'horizontalalignment', 'center', ...
                'backgroundcolor', xplr.colors('link_key', L.F.link_key));
            brick.controlpositions(L.h_label, L.h_p, [0, 1, 1, 0], [8, -21, -8-18, 18])
            % (close button)
            if ~new_figure
                x = brick.printnumber(ones(18), 'x', 'pos', 'center')';
                x(x == 1) = NaN;
                x = repmat(x, [1, 1, 3]);
                h_close = uicontrol('parent', L.h_p, 'cdata', x, 'callback', @(u,e)delete(L));
                brick.controlpositions(h_close, L.h_p, [1, 1], [-8-18, -3-18, 18, 18])
            end
            % (group button)
            ctrl = brick.propcontrol(L, 'sel_mult_in', 'togglebutton', ...
                {'parent', L.h_p, 'string', 'G'});
            brick.controlpositions(ctrl.hu, L.h_p, [0, 1], [8, -3-18, 18, 18])
            
            % context menu
            init_local_menu(L)
            
            % list and event (bottom-up)
            if brick.boolean(L.scroll_wheel)
                L.scroll_wheel = 'on'; % this will automaticall register scroll wheel
            end
            if isempty(get(L.h_f, 'WindowButtonMotionFcn'))
                % force update of current position when moving the mouse
                % around
                set(L.h_f, 'WindowButtonMotionFcn', @(u,e)do_nothing())
            end
            
            % watch filter
            function filter_changed(~, e)
                if strcmp(e.type, 'filter'), display_selection(L), end
            end
            brick.connect_listener(F, L, 'ChangedOperation', @filter_changed);
            
            % auto-delete
            set(L.h_list, 'deletefcn', @(u,e)delete(L))
            addlistener(F, 'ObjectBeingDestroyed', @(u,e)delete(L));

            % update display (here, just sets the correct value)
            preformat_values(L)
            display_selection(L)
            
            % set more properties
            if ~isempty(opt_add)
                set(L,opt_add{:})
            end
            
            % put object in base workspace for debugging purposes
            if nargin == 0, assignin('base', 'L', L), end
        end
        function init_local_menu(L)
            delete(L.menu)
            L.menu = uicontextmenu('parent', L.h_f);
            m = L.menu;
            set(L.h_list, 'UIContextMenu', m)
            
            uimenu(m, 'label', 'New singleton selections', 'callback', @(u,e)event(L, 'newuni'))
            uimenu(m, 'label', 'New group selection [A]', 'callback', @(u,e)event(L, 'newgroup'))
            uimenu(m, 'label', 'Add to selection', 'callback', @(u,e)event(L, 'add'))

            uimenu(m, 'label', 'Define new group...', 'callback', @(u,e)event(L, 'definegroup'), 'separator', 'on')

            uimenu(m, 'label', 'Sort selections according to list order', 'callback', @(u,e)event(L, 'sort'), 'separator', 'on')
            uimenu(m, 'label', 'Reorder selections...', 'callback', @(u,e)event(L, 'reorder'))
            
            uimenu(m, 'label', 'Remove highlighted group(s)', 'callback', @(u,e)event(L, 'rmgroup'), 'separator', 'on')
            uimenu(m, 'label', 'Remove all groups', 'callback', @(u,e)event(L, 'rmgroupall'))
            uimenu(m, 'label', 'Remove highlighted individuals', 'callback', @(u,e)event(L, 'rmuni'))
            uimenu(m, 'label', 'Remove all individuals', 'callback', @(u,e)event(L, 'rmuniall'))
            uimenu(m, 'label', 'Remove highlighted selections', 'callback', @(u,e)event(L, 'rm'))
            uimenu(m, 'label', 'Remove all selections', 'callback', @(u,e)event(L, 'rmall'))

            uimenu(m, 'label', 'Select all', 'separator', 'on', 'callback', @(u,e)event(L, 'selectall'))
            
            brick.propcontrol(L, 'sel_mult_in', ...
                {'menuval', {true, false}, {'individuals', 'group'}}, ...
                {'parent', m, 'label', 'Temporary selection', 'separator', 'on'});
            brick.propcontrol(L, 'selection_prompt_name', ...
                {'menu', {'all', 'groups', 'none'}, {'all selections', 'group selections only', 'none'}}, ...
                {'parent', m, 'label', 'Prompt for selection name'});
            
            % scroll wheel behavior: changing it would make sense only if
            % list is inside a figure with other elements, which does not
            % occur in XPLOR at the present time
            % and anyway, there are some bugs, and the whole
            % windowcallbackmanager thing needs to be replaced by calls to
            % iptaddcallback
            %             m1 = uimenu(m,'label','scroll wheel','separator','on');
            %             brick.propcontrol(L,'scroll_wheel', ...
            %                 {'menu', {'on' 'off' 'default'}}, ...
            %                 {'parent',m1,'label','Scroll wheel behavior'});
            %             uimenu(m1,'label','make default in figure', ...
            %                 'callback',@(u,e)set(L,'scroll_wheel','default'));

            % Load/save selections
            uimenu(m, 'label', 'Load selections...', 'separator', 'on', ...
                'callback', @(u,e)L.selection_load())
            uimenu(m, 'label', 'Save selections', 'enable', brick.onoff(~isempty(L.selection_save_file)), ...
                'callback', @(u,e)L.selection_save(L.selection_save_file))
            uimenu(m, 'label', 'Save selections as...', ...
                'callback', @(u,e)L.selection_save())
        end
        function selection_save(L, fname)
            if nargin<2 || isempty(fname)
                prompt = 'Select file for saving selections';
                fname = brick.savefile('*.xpls', prompt, L.selection_save_file);
                if isequal(fname,0), return, end
            end
            L.F.save_to_file(fname);
            L.selection_save_file = fname;
        end
        function selection_load(L, fname)
            if nargin<2
                fname = brick.getfile('*.xpls','Select selections file');
                if isequal(fname,0), return, end
            end
            try
                L.F.load_from_file(fname);
                L.selection_save_file = fname;
            catch ME
                errordlg(ME.message)
            end
        end
        function delete(L)
            delete@xplr.GraphNode(L)
            if ~isvalid(L) && ~isprop(L, 'h_list'), return, end
            brick.delete_valid(L.h_list, L.h_label)
        end
    end
       
    % Events
    methods
        function keypress(L, e)
            switch e.Key
                case 'a'
                    event(L, 'newgroup')
                case {'insert', 'delete'}
                    event(L, 'scroll', brick.switch_case(e.Key, 'insert', -1, 'delete', 1))
            end
        end
        function event(L, flag, varargin)
            % possible values for flag:
            % - select          selection by user click in list display
            % - unisel          selection by user double-click
            % - ('scroll',n)    scrolling
            % - newuni, newgroup
            % - definegroup
            % - add
            % - sort
            % - reorder...
            % - rmall, rmgroup, rmgroupall, rmuni, rmuniall, rm         
            
            % selected list entries
            val = get(L.h_list, 'value');
            %L.F.shared.list.cursel = val; % make available to all lists acting on this filter what is the new entries selection
            
            % get the current selections
            sel_inds = L.F.F.indices; % it is important here not to get L.F.indices, because we do not want the point selection to appear here
            n_sel = length(sel_inds);
            is_unis_el = false(1, n_sel);
            for i=1:n_sel, is_unis_el(i) = isscalar(sel_inds{i}); end % faster than calling brick.map
            
            % which selections are soft
            soft_sel = L.F.F.header_out.get_value('soft_selection'); % same as above
            soft_sel = [soft_sel{:}]; % faster than cell2mat
            if isempty(soft_sel), soft_sel = false(1,n_sel); end
            i_soft = find(soft_sel);
            n_soft = length(i_soft);
            i_solid = find(~soft_sel);
            n_solid = n_sel-n_soft;
            
            % modify flag if needed
            switch flag
                case 'select'
                    if strcmp(get(L.h_f, 'selectiontype'), 'open')
                        flag = 'unisel';
                    end
                case 'selectall'
                    set(L.h_list, 'value', 1:L.F.sz_in)
                    flag = 'select';
                case 'scroll'
                    n = varargin{1};
                    soft = sel_inds(soft_sel);
                    if length(soft) >= 2 && all(brick.map(@isscalar, soft)) ...
                            && all(ismember(diff([soft{:}]), [0, 1]))
                        % current soft selection consists of a range of
                        % consecutive values -> select the same number of
                        % values before or after
                        % Note that when hitting the lower or upper value,
                        % this value can be repeated several times in the
                        % selection in order to remember how many values
                        % where selected in total.
                        soft = [soft{:}];
                        n_val = length(soft);
                        if n < 0
                            % decreasing values
                            if any(diff(soft) == 0)
                                % if last value was selected several times
                                % show the n_val last values
                                val = soft(end) + (-n_val + 1:0);
                            else
                                % otherwise show the n_val values before
                                % soft(1)
                                val = soft(1) + n*n_val + (0:n_val-1);
                            end
                        else
                            % increasing values, same idea
                            if any(diff(soft) == 0)
                                val = soft(1) + (0:n_val - 1);
                            else
                                val = soft(end) + (n-1)*n_val + (1:n_val);
                            end
                        end
                        val = brick.coerce(val, 1, L.F.sz_in);
                    else
                        val = brick.coerce(L.F.index + n, 1, L.F.sz_in);
                        if val == L.F.index, return, end
                    end
                    set(L.h_list, 'value', val, 'listboxtop', val(1) - 2);
                    flag = 'select';
            end            
            
            % action (or only determine shich selections to remove and
            % which to add)
            idx_rm = [];
            new_sel = [];
            new_is_soft = false;
            switch flag
                case 'select'
                    if isscalar(val)
                        L.F.P.index_exact = val;
                    end
                    if isempty(val) || (isscalar(val) && isempty(i_solid)) ...
                            || (isscalar(val) && isequal({val}, sel_inds(i_soft)))
                        % remove temporary selection in the following
                        % cases:
                        % - user un_selected all list items
                        % - no solid selection and a single selected item
                        % - repeated selection of temporaray selection item
                        idx_rm = i_soft;
                    else
                        % new temporary selection
                        idx_rm = i_soft;
                        new_sel = build_current_selection(L, L.sel_mult_in);
                        new_is_soft = true;
                    end
                case 'unisel'
                    % double-click -> make new solid selection with current
                    % index, or remove it
                    if ~isscalar(val), return, end
                    k_uni_sel = find(is_unis_el);
                    idx_rm = k_uni_sel([sel_inds{k_uni_sel}] == val); % index of already-existing selection with this value
                    if isempty(idx_rm)
                        % create
                        new_sel = build_current_selection(L, true);
                    elseif soft_sel(idx_rm)
                        % make solid
                        new_sel = build_current_selection(L, true);
                    else
                        % remove
                    end
                case 'newuni'
                    idx_rm = i_soft;
                    new_sel = build_current_selection(L, true);
                case 'newgroup'
                    idx_rm = i_soft;
                    new_sel = build_current_selection(L, false);
                case 'definegroup'
                    str = inputdlg('Define selection', '', 1, {['1:', num2str(L.F.sz_in)]}); % TODO: continue!!!
                    if isempty(str), disp 'interrupted', return, end
                    try
                        val = evalin('base', ['[', str{1}, ']']);
                        if ~iscell(val), val = {val}; end
                        new_sel = xplr.SelectionND(length(val));
                        for i=1:length(val)
                            new_sel(i) = xplr.SelectionND(L.sel_type,val{i}, L.F.header_in.n);
                        end
                    catch
                        errordlg('Command could not be evaluated correctly')
                        return
                    end
                case 'add'
                    new_sel = xplr.SelectionND('point1D', val);
                    if n_soft, update_selection(L.F, 'remove', i_soft), end % remove all soft selections
                    if ~isempty(i_solid)
                        update_selection(L.F, 'add', i_solid(end), new_sel)
                    else
                        update_selection(L.F, 'new', new_sel)
                    end
                    return
                case 'sort'
                    % remove all soft selections
                    if n_soft, update_selection(L.F, 'remove', i_soft), end
                    sel_inds(i_soft) = [];
                    % reorder other selections
                    idx_first = brick.map(sel_inds, @(v)v(1), 'array');
                    [~, ord] = sort(idx_first);
                    update_selection(L.F, 'perm', ord)
                    return
                case 'reorder'
                    % first remove all soft selections
                    if n_soft, update_selection(L.F, 'remove', i_soft), end
                    n_sel = n_solid;
                    % prompt for reordering
                    ord = brick.input('new order', 1:n_sel);
                    if length(unique(ord)) < length(ord) || ~all(ismember(ord, 1:n_sel))
                        waitfor(errordlg('Not a valid permutation or subset'))
                        return
                    end
                    if length(ord) < n_sel
                        answer = questdlg('This is not a permutation: some selections will be removed','', ...
                            'OK','Cancel','OK');
                        if strcmp(answer, 'Cancel'), return, end
                        idx_rm = setdiff(1:n_sel, ord);
                        update_selection(L.F, 'remove', idx_rm)
                        [~, ord_sort] = sort(ord);
                        ord(ord_sort) = 1:length(ord);
                    end
                    update_selection(L.F, 'perm', ord)
                    return
                case 'rmall'
                    % set empty selections rather than remove all existing
                    % ones: this can performs some clean-up when errors
                    % occured previously
                    new_sel = xplr.SelectionND.empty(1, 0);
                    update_selection(L.F, 'all', new_sel)
                    return
                case {'rmgroup', 'rmgroupall', 'rmuni', 'rmuniall', 'rm'}
                    if strfind(flag, 'all')
                        range = 1:n_sel;
                        flag = strrep(flag, 'all', '');
                    else
                        range = brick.find(@(x)intersect(x, val), sel_inds);
                    end
                    rm_mask = false(1, n_sel);
                    switch flag
                        case 'rmgroup'
                            rm_mask(range) = ~is_unis_el(range);
                        case 'rmuni'
                            rm_mask(range) = is_unis_el(range);
                        case 'rm'
                            rm_mask(range) = true;
                    end
                    rm_mask(i_soft) = true; % in any case, remove all soft selections
                    idx_rm = find(rm_mask);
            end
            
            % prompt for name of new selections
            name_options = {};
            if ~strcmp(L.selection_prompt_name, 'none') && ~new_is_soft
                any_name = false;
                n_new = length(new_sel);
                names = cell(1, n_new);
                for i = 1:n_new
                    if strcmp(L.selection_prompt_name, 'groups') && isscalar(new_sel(i).data_ind), continue, end
                    name = inputdlg('Group name', 'xplor');
                    if isempty(name), continue, end
                    names{i} = name;
                    any_name = true;
                end
                if any_name
                    name_options = {'Name', names};
                end
            end
            
            % remove/change/add new selections
            n_new = length(new_sel);
            n_rm = length(idx_rm);
            if n_new == 0 && n_rm == 0
                % happens for example when there are no selection, and the
                % point is moved -> nothing to do
            elseif n_rm == 0
                update_selection(L.F, 'new', new_sel, 'soft_selection', new_is_soft)
            elseif n_new == 0
                update_selection(L.F, 'remove', idx_rm)
            elseif n_new == n_rm
                update_selection(L.F, 'chg', idx_rm, new_sel, 'soft_selection', new_is_soft, name_options{:})
            elseif n_new > n_rm
                update_selection(L.F, 'chg&new', {idx_rm, n_sel + (1:n_new - n_rm)}, new_sel, 'soft_selection', new_is_soft, name_options{:})
            elseif n_rm > n_new
                update_selection(L.F, 'chg&rm', {idx_rm(1:n_new), idx_rm(n_new + 1:n_rm)}, new_sel, 'soft_selection', new_is_soft, name_options{:})
            end
        end
        function sel = build_current_selection(L, do_mult_in)
            val = get(L.h_list, 'value');
            if isempty(val), sel = []; return, end
            if do_mult_in
                n_sel = length(val);
                sel = xplr.SelectionND(n_sel);
                for i=1:n_sel
                    sel(i) = xplr.SelectionND(L.sel_type, val(i), L.F.header_in.n);
                end
            elseif L.F.header_in.is_measure && ~isscalar(val) && all(diff(val) == 1)
                % selection is a segment rather than a mere list of points
                sel = xplr.SelectionND('line1D', val([1, end]) + [-.5, .5], L.F.header_in.n);
            else
                sel = xplr.SelectionND(L.sel_type, val, L.F.header_in.n);
            end
        end
    end
        
    % Get/Set - scroll wheel
    methods
        function set.scroll_wheel(L, flag)
            switch flag
                case 'on'
                    brick.scroll_wheel_register(L.h_list, @(n)event(L, 'scroll', n)) %#ok<MCSUP>
                    L.scroll_wheel = 'on';
                case 'default'
                    brick.scroll_wheel_register(L.h_list, @(n)event(L, 'scroll', n), 'default') %#ok<MCSUP>
                    L.scroll_wheel = 'on';
                case 'off'
                    brick.scroll_wheel_register(L.h_list, flag) %#ok<MCSUP>
                    L.scroll_wheel = 'off';
                otherwise
                    error 'scroll_wheel value must be ''off'', ''on'' or ''default'''
            end
        end
    end
    
    % Get/Set
    methods
        function set.sel_mult_in(L, val)
            if val == L.sel_mult_in, return, end
            L.sel_mult_in = val;
            % update selection
            event(L, 'select')
        end
    end
    
    % Display
    methods (Access = 'private')
        function preformat_values(L)
            L.value_str = L.F.header_in.get_item_names();
        end
        function display_selection(L)
            % init list with names of items
            str = L.value_str;
            
            sel_inds = L.F.F.indices; % it is important here not to get L.F.indices, because we do not want the point selection to appear here
            n_sel = length(sel_inds);
            is_unis_el = false(1, n_sel);
            for i=1:n_sel, is_unis_el(i) = isscalar(sel_inds{i}); end % faster than calling brick.map
            soft_sel = L.F.F.header_out.get_value('soft_selection'); % same as above
            soft_sel = [soft_sel{:}]; % faster than cell2mat
            if isempty(soft_sel), soft_sel = false(1, n_sel); end
                        
            % mark selections
            for k_sel = find(is_unis_el & ~soft_sel)
                ind = sel_inds{k_sel};
                str{ind} = [str{ind}, '[', num2str(k_sel), ']'];
            end
            for k_sel = find(~is_unis_el & ~soft_sel)
                for ind = sel_inds{k_sel}
                    str{ind} = [str{ind}, '[group', num2str(k_sel), ']'];
                end
            end
            for k_sel = find(is_unis_el & soft_sel)
                ind = sel_inds{k_sel};
                str{ind} = [str{ind}, '*'];
            end
            for k_sel = find(~is_unis_el & soft_sel)
                for ind = sel_inds{k_sel}
                    str{ind} = [str{ind}, '*g'];
                end
            end
            
            % update display!
            top = get(L.h_list, 'listboxtop'); % the portion of the list that is shown moves when setting the string -> reset it to its current position
            set(L.h_list, 'string', str, 'ListboxTop', top)
            if n_sel == 0
                % if no selection, highlight the point selection
                set(L.h_list, 'value', L.F.point_index)
            else
                set(L.h_list, 'value', [sel_inds{soft_sel}])
            %             elseif isfield(L.F.shared,'list')
            %                 set(L.h_list,'value',L.F.shared.list.cursel)
            end
        end
            
    end
    
end
     
%---
function do_nothing()
% setting this function as main figure WindowButtonMotionFcn forces the
% figure to update CurrentPoint whenever the mouse is moved, and therefore
% have this property set currently even when clicking active controls
% (unfortunately, clicking a control does not set this property)
end
