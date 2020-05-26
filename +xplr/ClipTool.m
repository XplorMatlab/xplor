classdef ClipTool < xplr.GraphNode
% ClipTool 
    
    % properties controlled by (and as seen from) the menu
    properties (SetObservable=true, AbortSet, SetAccess='private')
        auto_clip_mode = 'minmax';
    end
    properties (SetObservable=true, AbortSet)
        auto_clip_mode_no_center = 'minmax'     % 
        center = []                         % [], 0 or 1
        adjust = 'none'                     % 'none', 'mean'
        span = 'local'                      % 'curview', 'local' or 'shared (share_name)'
    end

    properties (SetAccess='private')
        share_name = ''                      % 
        menu                                % 
    end
    
    events
        ChangedClip
    end
   
    methods
        function C = ClipTool(hf)
            % constructor
            if nargin==0
                hf = figure(826);
                clf(hf)
                set(hf,'handlevisibility', 'off')
            end
            build_menu(C, hf)
        end
        function delete(C)
            delete@xplr.GraphNode(C)
            if ~isprop(C, 'menu'), return, end
            delete_valid(C.menu)
        end
        function build_menu(C,hf)
            if isempty(C.menu)
                m = uimenu('parent', hf, 'label', 'Clipping');
                C.menu = m;
            else
                m = C.menu;
                if nargin>=2 && get(m, 'parent') ~= hf, error 'existing menu does not have the same parent figure', end
                delete(get(m, 'children'))
            end
            % auto-clip specifications
            % (create a submenu whose label will update automatically)
            P = fn_propcontrol(C, 'auto_clip_mode', ...
                {'menuval', {}}, ...
                m, 'label', 'Auto-Clip Method');
            m1 = P.hparent;
            % (some possible values)
            fn_propcontrol(C, 'auto_clip_mode_no_center', ...
                {'menugroup', {'minmax', 'std1', 'std2', 'std3', 'std5', 'std.5', 'std.2', 'prc.1', 'prc1', 'prc5'}, ...
                {'Min to Max', '1 STD', '2 STD', '3 STD', '5 STD', '1/2 STD', '1/5 STD', '[.1% 99.9%]', '[1% 99%]', '[5% 95%]'}}, ...
                m1);
            % (all possibilities through additional input dialogs)
            uimenu(m1, 'label', 'Use Mean and STD...', 'separator','on', 'callback', @(u,e)set_auto_clip_mode(C, 'setstd'))
            uimenu(m1, 'label', 'Use Percentiles...', 'callback', @(u,e)set_auto_clip_mode(C, 'setprc'))
            % (fix the center)
            fn_propcontrol(C, 'center', ...
                {'menugroup', {0, 1, []}, {'center on 0', 'center on 1'}}, ...
                m1);
            % set clip
            uimenu(m, 'label', 'Set Clip Manually...', 'callback', @(u,e)set_manual_clip(C))
            uimenu(m, 'label', 'Do Auto-Clip', 'callback', @(u,e)notify(C, 'ChangedClip', xplr.EventInfo('clip', 'automode')))
            % adjust
            fn_propcontrol(C, 'adjust', ...
                {'menuval', {'none', 'mean'}, {'None', 'Adjust each Element by its Mean Value'}, {'None', 'Mean'}}, ...
                m, 'label', 'Adjust', 'separator', 'on');
            % span
            fn_propcontrol(C, 'span', ...
                {'menuval', {'curview', 'local', 'chg_share'}, {'Local, Adjust to Current View', 'Local', 'Shared...'}, {'Current View', 'Local', ''}}, ...
                m, 'label', 'Span');
        end
    end
    
    % Setting clip manually
    methods
        function set_manual_clip(C, clip)
            if nargin < 2
                clip = fn_input('Enter min and max', [0, 1]);
                if length(clip) ~= 2, return, end % value not valid, do not use it
            end
            notify(C, 'ChangedClip', xplr.EventInfo('clip', 'clip', clip))
        end
    end
    
    % Setting autoclip mode
    methods
        function set_auto_clip_mode(C, comp)
            switch comp
                case 'setstd'
                    nstd = fn_input('Mean +/- N*std', 1, 'stepper 1 0 Inf');
                    if isempty(nstd) || nstd == 0, return, end
                    C.auto_clip_mode_no_center = ['std', num2str(nstd)];
                case 'setprc'
                    prc = fn_input('Percentiles (one or two values between 0 and 100)', 1);
                    if isempty(prc), return, end
                    switch length(prc)
                        case 1
                            C.auto_clip_mode_no_center = ['prc', num2str(prc)];
                        case 2
                            C.auto_clip_mode_no_center = ['prc', num2str(prc(1)), '-', num2str(prc(2))];
                        otherwise
                            % value not valid, do not set property
                    end
                otherwise
                    C.auto_clip_mode_no_center = comp;
                    return
            end
        end
        function set.auto_clip_mode_no_center(C, comp)
            % set property
            C.auto_clip_mode_no_center = comp;
            % add the centering information to final autoclip mode
            C.autoclipmode = C.auto_clip_mode_no_center;
            if ~isempty(C.center), C.autoclipmode = [C.autoclipmode, '[', num2str(C.center), ']']; end
            % notify the change
            notify(C, 'ChangedClip', xplr.EventInfo('clip', 'automode'))
        end
        function set.center(C, center)
            % set property
            C.center = center;
            % update autoclip mode
            C.auto_clip_mode = C.auto_clip_mode_no_center;
            if ~isempty(C.center), C.auto_clip_mode = [C.auto_clip_mode, '[', num2str(C.center), ']']; end
            % notify the change
            notify(C, 'ChangedClip', xplr.EventInfo('clip', 'automode'))
        end
    end
    
    % Adjusting
    methods
        function set.adjust(C,adj)
            C.adjust = adj;
            notify(C, 'ChangedClip', xplr.EventInfo('clip', 'adjust'))
        end
    end
    
    % Sharing
    methods
        function set.span(C,span)
            if strcmp(span, 'chg_share')
                C.share_name = fn_input('Clip sharing name', '#1');
                C.span = ['Shared (', C.share_name, ')']; %#ok<*MCSUP>
            else
                C.span = span;
            end
            notify(C, 'ChangedClip', xplr.EventInfo('clip', 'span'))
            % shared clipping not really handled yet
        end
    end
end
