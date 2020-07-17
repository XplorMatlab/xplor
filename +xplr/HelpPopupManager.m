classdef HelpPopupManager < matlab.mixin.SetGet
    
    % TODO:
    % - complete methods popup_window and button_clicked
    % - add a 'Help' menu in View.m with items "Show help popups" (use
    %   function fn_prop_control on the HelpPopupManager object to handle
    %   the check mark) and "Reset popup display" 
    
    properties
        % "global" properties
        displayed_identifiers_session = {};
        displayed_identifiers_disk = {};
        
        % properties related to the currently opened popup
        current_identifier
        current_popup           % probably an uifigure
        current_object
        current_object_state
    end
    properties (SetObservable)
        do_show_popups = false
    end
    
    % Help popup manager
    methods (Access='private')
        % Constructor
        function self = HelpPopupManager()
            % read from disk properties that are persistent accross Matlab
            % sessions
            self.load_from_disk('do_show_popups')
            self.load_from_disk('displayed_identifiers_disk')
        end
        function load_from_disk(self, prop)
            % loadprop
            fsave = fn_userconfig('configfolder', 'xplor_popup');
            warning('off', 'MATLAB:load:variableNotFound')
            try %#ok<TRYNC>
                self.(prop) = fn_loadvar(fsave, prop);
            end
            warning('on', 'MATLAB:load:variableNotFound')
        end
        function save_to_disk(self, prop)
            fsave = fn_userconfig('configfolder', 'xplor_popup');
            s = struct(prop, {self.(prop)});
            if exist(fsave, 'file')
                save(fsave, '-STRUCT', 's', '-APPEND');
            else
                save(fsave, '-STRUCT', 's');
            end
        end
    end
    methods (Static)
        function manager = get_popup_manager()
            % Unique popup manager is attached to the root graphic object.
            % This is preferrable to using a global variable that might be
            % deleted with the 'clear' command.
            M0 = getappdata(0, 'xplor_popup');
            if isempty(M0)
                M0 = xplr.HelpPopupManager();
                setappdata(0, 'xplor_popup', M0)
            end
            manager = M0;
        end
        function build_menu(hf)
%             manager = xplr.HelpPopupManager.get_popup_manager();% test fn_prop
%             fn_propcontrol(manager,'do_show_popups','menu',...  % test fn_prop
%                {'parent',hf,'label','Help'});                   % test fn_prop
            menu = uimenu('parent',hf,'label','Help');
            uimenu(menu,'label','Show help popups');
            uimenu(menu,'label','Reset popups display');
        end
    end
    
    % Identifiers lists
    methods (Access='private')
        function value = was_identifier_displayed(self, identifier)
            value = ismember(identifier, self.displayed_identifiers_disk) ...
                || ismember(identifier, self.displayed_identifiers_session);
        end
        function add_to_disk_list(self, identifier)
            self.displayed_identifiers_disk{end+1} = identifier;
            self.save_to_disk('displayed_identifiers_disk');
        end
        function add_to_session_list(self, identifier)
            self.displayed_identifiers_session{end+1} = identifier;
        end
        function reset_identifier_lists(self)
            self.displayed_identifiers_session = {};
            self.displayed_identifiers_disk = {};
            self.save_to_disk('displayed_identifiers_disk');
        end
    end
    
    % Do show popup
    methods (Static)
        function value = get_do_show_popups()
            manager = xplr.HelpPopupManager.get_popup_manager();
            value = manager.do_show_popups;
        end
        function set_do_show_popups(value)
            manager = xplr.HelpPopupManager.get_popup_manager();
            manager.do_show_popups = value;
        end
    end
    
    % Popup window
    methods (Static)
        function popup_window(file_name, object)
            % create dialog window with appropriate buttons
            % - display message = web page
            % - highlight a control (stop highlight when window closed)
            % - buttons: "Do not show again" (=default), "Show me again later"
            
            manager = xplr.HelpPopupManager.get_popup_manager();

            % For testing purposes
            manager.reset_identifier_lists()
            disp reset
            
            ~manager.do_show_popups
            manager.was_identifier_displayed(file_name)
            % is this message in a list of already displayed messages?
            if ~manager.do_show_popups || manager.was_identifier_displayed(file_name)
                return
            end

            if nargin >= 2, manager.object = object; end
            
            manager.current_identifier = file_name;
            
            % for the moment, just display file_name
            % click
            disp(file_name)
            
            % Test if Matlab version is earlier than R2019b
            uihtml_ok = ~verLessThan('matlab','9.7');
            
            if uihtml_ok
                % Create figure
                popup_figure = uifigure;
                popup_figure.Position = [500 500 380 445];
                % Add HTML to figure
                message = uihtml(popup_figure);
                message.Position = [10 10 360 420];
                message.HTMLSource = file_name;
                % Add buttons
                again_button = uibutton(popup_figure);
                again_button.Text = 'Show me again later';
                again_button.Position = [25 50 150 30];
                not_again_button = uibutton(popup_figure);
                not_again_button.Text = 'Don''t show me again';
                not_again_button.Position = [200 50 150 30];
            else
                disp('Matlab version earlier than R2019b detected')
                % Read HTML file and delete tags to keep raw text only
                popup_html_text = fileread(file_name);
                popup_text = regexprep(popup_html_text, '<.*?>', '');
                % Questdialog popup with raw text and 2 buttons
                questdlg(popup_text, 'Help box', 'Show me again later', 'Don''t show me again later', 'Show me again later')
            end
            
            % Emulate button click
            manager.button_clicked(true);
        end
    end
    methods
        function button_clicked(self, do_not_show_again)
            
            % close window, restore object previous state
            
            % store identifier
            if do_not_show_again
                self.add_to_disk_list(self.current_identifier)
            else
                self.add_to_session_list(self.current_identifier)
            end
            self.current_identifier = [];
        end
    end

end
