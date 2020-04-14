classdef HelpPopupManager < handle
    
    % TODO:
    % - complete methods popup_window and button_clicked
    % - add a 'Help' menu in view.m with items "Show help popups" (use
    %   function fn_prop_control on the HelpPopupManager object to handle
    %   the check mark) and "Reset popup display" 
    
    properties
        % "global" properties
        do_show_popups = false;
        displayed_identifiers_session = {};
        displayed_identifiers_disk = {};
        
        % properties related to the currently opened popup
        current_identifier
        current_popup           % probably an uifigure
        current_object
        current_object_state
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
        function load_from_disk(self,prop)
            % loadprop
            fsave = fn_userconfig('configfolder','xplor_popup');
            warning('off','MATLAB:load:variableNotFound')
            try %#ok<TRYNC>
                self.(prop) = fn_loadvar(fsave,prop);
            end
            warning('on','MATLAB:load:variableNotFound')
        end
        function save_to_disk(self,prop)
            fsave = fn_userconfig('configfolder','xplor_popup');
            s = struct(prop,{self.(prop)});
            if exist(fsave,'file')
                save(fsave,'-STRUCT','s','-APPEND');
            else
                save(fsave,'-STRUCT','s');
            end
        end
    end
    methods (Static)
        function manager = get_popup_manager()
            % Unique popup manager is attached to the root graphic object.
            % This is preferrable to using a global variable that might be
            % deleted with the 'clear' command.
            M0 = getappdata(0,'xplor_popup');
            if isempty(M0)
                M0 = xplr.HelpPopupManager();
                setappdata(0,'xplor_popup',M0)
            end
            manager = M0;
        end
    end
    
    % Identifiers lists
    methods (Access='private')
        function value = was_identifier_displayed(self, identifier)
            value = ismember(identifier,self.displayed_identifiers_disk) ...
                || ismember(identifier,self.displayed_identifiers_session);
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
        function popup_window(filename, object)
            % create dialog window with appropriate buttons
            % - display message = web page
            % - highlight a control (stop highlight when window closed)
            % - buttons: "Do not show again" (=default), "Show me again later"
            
            manager = xplr.HelpPopupManager.get_popup_manager();
            
            % is this message in a list of already displayed messages?
            if ~manager.do_show_popups || manager.was_identifier_displayed(filename)
                return
            end

            if nargin>=2, manager.object = object; end
            
            manager.current_identifier = filename;
            
            % for the moment, just display filename and emulate button
            % click
            disp(filename)
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