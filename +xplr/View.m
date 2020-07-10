classdef View < xplr.GraphNode
    % Main window
    %
    
    % Components 
    properties (SetAccess = 'private')
        hf              % figure
        slicer          % operation to get 'slice' from 'data'
        C               % control of data operation
        D               % main display
        panels          % panels for: display, control, lists, + some buttons
        menu_xplor
        menu_help
    end
    properties (SetObservable = true)
        control_visible = false;  % logical - are the controls visible
    end
    
    % Dependent properties
    properties (Dependent, SetAccess = 'private')
        data
        slice
    end
        
    % Constructor
    methods
        function V = View(data, varargin)
            % function V = view(data[,optional slicer arguments...])

            % define headers if needed
            if ~isa(data, 'xplr.XData')
                head = xplr.EditHeader(data);
                if isempty(head)
                    % header editing was canceled
                    % (make Dependent properties 'data' and 'slice' valid
                    % to avoid errors when V will be displayed)
                    V.slicer = struct('data', [], 'slice', []);
                    return
                end
                data = xplr.XData(data, head);
            end
            
            % SLICER
            if ~isa(data,'xplr.XData'), error 'data argument must be a xplr.xdata object', end
            V.slicer = V.add_component(xplr.Slicer(V, data, varargin{:}));
            
            % LINKS
            % register view object to the bank
            xplr.Bank.register_view(V)

            % PANELS
            % open figure and create panels
            init_panels(V)
            
            % MENU XPLOR
            V.menu_xplor = uimenu('parent', V.hf, 'label', 'XPLOR', ...
                'callback', @(u,e)V.xplor_menu());
            
            % DISPLAY
            % (note: calling add_component will cause V.D to be deleted when
            % V will be deleted and vice-versa; the output of add_component
            % is simply its input)
            V.D = V.add_component(xplr.ViewDisplay(V));
            
            % CONTROL
            V.C = V.add_component(xplr.ViewControl(V));
            
            % MENU HELP
            V.menu_help = uimenu('parent', V.hf, 'label', 'Help', ...
                'callback', @(u,e)V.help_menu());
            
            % save object in base workspace
            assignin('base', 'V', V)
        end
        function delete(V)
            if ~isprop(V,'hf'), return, end
            delete_valid(V.hf, V.C, V.D)
        end
    end
    
    % Reset display
    methods
        function reset_display(V)
            delete(allchild(V.panels.display))
            V.D = xplr.ViewDisplay(V);
        end
    end
    
    % Panels
    methods
        function init_panels(V)
            % figure
            V.hf = figure('integerhandle', 'off', 'handlevisibility', 'off', 'visible', 'off', ...
                'numbertitle', 'off', 'name', 'XPLOR', 'tag', 'XPLOR', ...
                'menubar', 'none', ...
                'color', 'white');
            delete(findall(V.hf, 'parent', V.hf))
            set(V.hf, 'DeleteFcn', @(u,e)delete(V))
            set(V.hf, 'WindowButtonMotionFcn', @(u,e)do_nothing())
            if ~isempty(V.data.name), set(V.hf, 'name', ['XPLOR: ', V.data.name]), end
            
            % by defaults, callbacks should not be interruptible
            set(V.hf,'Interruptible', 'off', ...
                'DefaultUicontrolInterruptible', 'off', ...
                'DefaultLineInterruptible', 'off', ...
                'DefaultTextInterruptible', 'off')
            
            % panels: use panelorganizer tool
            V.panels.main = panelorganizer(V.hf, 'H', 2, [0, 1],[0, 1]); % left part has fixed width, right part adjusts automatically
            set(V.panels.main, 'bordermode', 'push')
            V.panels.display = V.panels.main.setSubPanel(2);
            V.panels.all_controls = V.panels.main.setSubOrg(1, 'V', 2,[1, 1], [1, 0]);
            connect_listener(V.panels.all_controls.hobj, V, 'Visible', 'PostSet', @(u,e)set(V, 'control_visible', get(V.panels.all_controls.hobj, 'visible')));
            V.panels.controls_width = 120;
            V.panels.control = V.panels.all_controls.setSubPanel(1);
            V.panels.list_combo = V.panels.all_controls.setSubOrg(2, 'H');
            V.control_visible = true; % start with controls visible; it is necessary to explicitely set property to true (rather than having true as default value) to update display correctly
            
            % switch of control visibility
            ppi = get(0, 'screenpixelsPerInch');
            control_on = fn_propcontrol(V, 'control_visible', {'pushbutton', {true, false}, {'>> Controls', '<< Controls'}}, ...
            	'parent', V.panels.display, 'backgroundcolor', 'w', 'fontsize', 800/ppi);
            fn_controlpositions(control_on.hu, V.panels.display, [0, 1], [2, -17, 65, 15])
            
            % now figure can be visible
            set(V.hf, 'visible', 'on')
        end
        function set.control_visible(V, b)
            if strcmp(b, 'toggle')
                b = ~V.control_visible;
            else
                b = boolean(b);
                if b == V.control_visible, return, end
            end
            % set property
            V.control_visible = b;
            % positions
            main = V.panels.main;
            if b && main.extents(1) == 0
                main.pushExtent(1, V.panels.controls_width, 'figleft'); %#ok<*MCSUP>
            elseif ~b && main.extents(1) > 0
                V.panels.controls_width = main.extents(1); % memorize for next controls re-opening
                main.pushExtent(1, 0, 'figleft');
            end
        end
        function show_controls(V)
            V.control_visible = true;
        end
        function hide_controls(V)
            V.control_visible = false;
        end
    end
    
    % Menu
    methods
        function xplor_menu(V)
            m = V.menu_xplor;
            delete(get(m,'children'))
            
            % Access View from command line
            uimenu(m, 'label', 'Access window from command line...', ...
                'callback', @(u,e)V.command_line_access())
            
            % Close XPLOR window(s)
            fig_name = get(V.hf, 'name');
            uimenu(m, 'label', ['Close all ''' fig_name ''' windows'], 'separator', 'on', ...
                'callback', @(u,e)close(findall(0, 'type', 'figure', 'name', fig_name)))
            uimenu(m, 'label', 'Close all XPLOR windows', ...
                'callback', @(u,e)close(findall(0, 'type', 'figure', 'tag', 'XPLOR')))
        end
        function command_line_access(V)
            persistent do_ask_name
            % suggest a variable name
            idx = 0;
            ok = false;
            while ~ok
                idx = idx + 1;
                name = ['V' num2str(idx)];
                if ~evalin('base',['exist(''' name ''', ''var'')'])
                    ok = true;
                else
                    V_name = evalin('base', name);
                    ok = (~isvalid(V_name) || V == V_name);
                end    
            end
            % prompt user for variable name
            if isempty(do_ask_name), do_ask_name = true; end
            if do_ask_name
                s = fn_structedit( ...
                    struct('Variable__name', name, 'do__not__ask__again', false));
                if isempty(s), return, end
                name = s.Variable__name;
                do_ask_name = ~s.do__not__ask__again;
            end
            % save 
            assignin('base', name, V)
            disp(['View object has been stored in variable ''' name '''.'])
        end
        function help_menu(V)
            m = V.menu_help;
            delete(get(m,'children'))
            
            popup_manager = xplr.HelpPopupManager.get_popup_manager();
            
            % Create Help menu
            uimenu(m, 'label', 'TEST', ...
                'callback', @(u,e)V.test_help())
            
            % Button Show help popups
            fn_propcontrol(popup_manager, 'do_show_popups', 'menu', ...
                {'parent', m, 'label', 'Show help popups'});
            
            % Button Reset popups display
            uimenu(m, 'label', 'Reset display', 'separator', 'on', ...
                'callback', @(u,e)popup_manager.reset_identifier_lists())
            
        end
        function test_help(V)
            disp("TEST HELP")
            disp(xplr.HelpPopupManager.get_popup_manager())
        end
    end
    
    % Dependent properties
    methods
        function data = get.data(V)
            data = V.slicer.data;
        end
        function slice = get.slice(V)
            slice = V.slicer.slice;
        end
    end
    
    % Handling of the axis
    methods (Access='private')
        function init_figure(V,in)
            if isempty(in), in = gca; end
            if strcmp(get(in, 'type'), 'axes')
                V.ha = in;
                V.hf = get(V.ha, 'parent');
            else
                fn_isfigurehandle(in) % create figure if in is a positive integer but not yet a figure handle
                V.hf = in;
                V.ha = axes('parent', V.hf);
            end
        end
        function init_axes(V)
            % structure with all secondary graphic objects
            s = struct;
            
            % axes on the side serve to generate "outside" mouse events
            col = get(V.hf, 'color')*.9;
            s.ha2 = axes('parent', V.hf, 'buttondownfcn', @(u,e)Mouse(V, 'outsidedown'));
            fn_controlpositions(s.ha2(1), V.ha, [0, 0, 1, 0], [0, -20, 0, 20]);
            s.ha2(2) = axes('parent', V.hf, 'buttondownfcn', @(u,e)Mouse(V, 'outsideleft'));
            fn_controlpositions(s.ha2(2), V.ha, [0, 0, 0, 1], [-20, 0, 20, 0]);
            s.ha2(3) = axes('parent', V.hf, 'buttondownfcn', @(u,e)Mouse(V, 'outsideboth'));
            fn_controlpositions(s.ha2(3), V.ha, [0, 0, 0, 0], [-20, -20, 20, 20]);
            set(s.ha2, 'handlevisibility', 'off', 'color', col, ...
                'xtick', [], 'xcolor', col, 'ytick', [], 'ycolor', col)
            uistack(s.ha2, 'bottom')
            
            % store structure
            V.hobj = s;
        end
        function init_pointer(V)
            h(1) = line('Parent', V.ha, 'xdata', [0, 0]);
            h(2) = line('Parent', V.ha, 'ydata', [0, 0]);
            h(3) = line('Parent', V.ha, 'xdata', 0, 'ydata', 0, 'marker', '.', 'linestyle', 'none'); % a single point
            V.hobj.pointer = h;
            if V.isimage
                set(h, 'Color', 'white')
            else
                set(h, 'Color', 'k')
                set(h(2:3), 'visible', 'off')
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
