classdef ClipTool < xplr.GraphNode
% ClipTool
    
    % properties controlled by (and as seen from) the menu
    properties (SetAccess='private')
        D
    end
    properties (SetObservable=true, AbortSet, SetAccess='private')
        auto_clip_mode = 'minmax'
    end
    properties (SetObservable=true, AbortSet)
        auto_clip_mode_no_center = 'minmax'     %
        center = []                         % [], 0 or 1
        independent_dim_id_mem = []          % dimension ID of dimensions along which clipping is not uniform
        align_signals = ''                  % '', 'nmean', 'nmedian'
        adjust_to_view = true
        buttons_only_for_current_cells = false
    end
    properties (Dependent, SetAccess='private')
        independent_dim
        linked_dim
    end

    properties (SetAccess='private')
        share_name = ''                      %
        menu                                % 
    end
    
    % Constructor, menu
    methods
        function C = ClipTool(D)
            % constructor
            C.D = D;
            C.menu = uimenu('parent', D.V.hf, 'label', 'Clipping', ...
                'callback', @(u,e)clip_menu(C));
        end
        function delete(C)
            delete@xplr.GraphNode(C)
            if ~isprop(C, 'menu'), return, end
            brick.delete_valid(C.menu)
        end
        function clip_menu(C)
            m = C.menu;
            delete(get(m,'children'))
            
            % auto-clip specifications
            % (create a submenu whose label will update automatically)
            P = brick.propcontrol(C, 'auto_clip_mode', ...
                {'menuval', {}}, ...
                m, 'label', 'Auto-Clip Method');
            m1 = P.hparent;
            % (some possible values)
            brick.propcontrol(C, 'auto_clip_mode_no_center', ...
                {'menugroup', {'minmax', 'std1', 'std2', 'std3', 'std5', 'std.5', 'std.2', 'prc.1', 'prc1', 'prc5'}, ...
                {'Min to Max', '1 STD', '2 STD', '3 STD', '5 STD', '1/2 STD', '1/5 STD', '[.1% 99.9%]', '[1% 99%]', '[5% 95%]'}}, ...
                m1);
            % (all possibilities through additional input dialogs)
            uimenu(m1, 'label', 'Use Mean and STD...', 'separator','on', 'callback', @(u,e)set_auto_clip_mode(C, 'setstd'))
            uimenu(m1, 'label', 'Use Percentiles...', 'callback', @(u,e)set_auto_clip_mode(C, 'setprc'))
            % (fix the center)
            brick.propcontrol(C, 'center', ...
                {'menugroup', {0, 1, []}, {'center on 0', 'center on 1'}}, ...
                m1);
            
            % auto-clip
            if ~isempty(C.independent_dim)
                uimenu(m, 'label', 'Do Auto-Clip (current cell(s) only)', 'callback', @(u,e)C.D.auto_clip(false))
                uimenu(m, 'label', 'Do Auto-Clip (all cells)', 'callback', @(u,e)C.D.auto_clip(true))
                brick.propcontrol(C, 'buttons_only_for_current_cells', 'menu', ...
                    {'parent', m, 'label', 'Use clip buttons to control only current cell(s)'});
            else
                uimenu(m,'label', 'Do Auto-Clip', 'callback', @(u,e)C.D.auto_clip(true))
            end
            
            % adjust options
            brick.propcontrol(C, 'adjust_to_view', 'menu', ...
                {'parent', m, 'label', 'Always adjust clipping to current view', 'separator', 'on'});
            clip_dim = C.D.clip_dim;
            if ~isempty(clip_dim)
                if isscalar(clip_dim)
                    checked = ismember(clip_dim, C.independent_dim);
                    uimenu(m, 'label', 'Independent clipping range for each grid cell', ...
                        'checked', checked, ...
                        'callback', @(u,e)set_independent_dim(C, clip_dim, ~checked))
                else
                    m1 = uimenu('parent', m, 'label', 'Independent clipping range for dimension(s)');
                    for dimID = clip_dim
                        checked = ismember(dimID, C.independent_dim);
                        uimenu(m1, 'label', C.D.slice.dimension_label(dimID), ...
                            'checked', checked, ...
                            'callback', @(u,e)set_independent_dim(C, dimID, ~checked))
                    end
                    uimenu(m1, 'label', '(all)', 'separator', 'on', ...
                        'callback', @(u,e)set_independent_dim(C, clip_dim, true))
                    uimenu(m1, 'label', '(none)', ...
                        'callback', @(u,e)set_independent_dim(C, clip_dim, false))
                end
            end
            if strcmp(C.D.display_mode, 'time courses')
                brick.propcontrol(C, 'align_signals', ...
                    {'menuval', {'', 'nmean', 'nmedian'}, ...
                    {'(do not center)', 'on their mean', 'on their median'}, {'(none)', 'mean', 'median'}}, ...
                    {'parent', m, 'label', 'Align signals'});
            end
            % manual clip
            if ~isempty(C.independent_dim)
                uimenu(m, 'label', 'Set Clip Manually... (current cell(s) only)', 'separator', 'on', ...
                    'callback', @(u,e)C.set_manual_clip([], false))
                uimenu(m, 'label', 'Set Clip Manually... (all cells)', ...
                    'callback', @(u,e)C.set_manual_clip([], true))
            else
                uimenu(m, 'label', 'Set Clip Manually...', 'separator', 'on', ...
                    'callback', @(u,e)C.set_manual_clip([], true))
            end
        end
    end

    % Get dependent
    methods
        function dim = get.independent_dim(C)
            indp_dim = C.D.slice.dimension_number(C.independent_dim_id_mem);
            dim = intersect(C.D.clip_dim, indp_dim);
        end
        function dim = get.linked_dim(C)
            indp_dim = C.D.slice.dimension_number(C.independent_dim_id_mem);
            dim = setdiff(C.D.clip_dim, indp_dim);
        end
    end
    
    % Setting clip manually
    methods
        function set_manual_clip(C, clip, all_cells)
            if nargin < 2 || isempty(clip)
                clip = brick.input('Enter min and max', [0, 1]);
                if length(clip) ~= 2, return, end % value not valid, do not use it
            end
            % update display
            if nargin<3, all_cells = true; end
            C.D.set_clip(clip, all_cells)
        end
    end
    
    % Setting autoclip mode
    methods
        function set_auto_clip_mode(C, comp)
            switch comp
                case 'setstd'
                    nstd = brick.input('Mean +/- N*std', 1, 'stepper 1 0 Inf');
                    if isempty(nstd) || nstd == 0, return, end
                    C.auto_clip_mode_no_center = ['std', num2str(nstd)];
                case 'setprc'
                    prc = brick.input('Percentiles (one or two values between 0 and 100)', 1);
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
            C.auto_clip_mode = C.auto_clip_mode_no_center;
            if ~isempty(C.center), C.auto_clip_mode = [C.autoclipmode, '[', num2str(C.center), ']']; end
            % update display
            C.D.auto_clip(true)
        end
        function set.center(C, center)
            % set property
            C.center = center;
            % update autoclip mode
            C.auto_clip_mode = C.auto_clip_mode_no_center;
            if ~isempty(C.center), C.auto_clip_mode = [C.auto_clip_mode, '[', num2str(C.center), ']']; end
            % update display
            C.D.auto_clip(true)
        end
    end
    
    % Adjusting
    methods
        function set_independent_dim(C, dim_id, independent)
            if nargin<3, independent = true; end
            dim_id = C.D.slice.dimension_id(dim_id);
            if independent
                C.independent_dim_id_mem = union(C.independent_dim_id_mem, dim_id);
            else
                C.independent_dim_id_mem = setdiff(C.independent_dim_id_mem, dim_id);
            end
            % update clip if some dimensions are no more independent
            C.D.auto_clip(true)
        end
        function set.adjust_to_view(C, value)
            C.adjust_to_view = value;
            % update display
            if value
                C.D.auto_clip(true)
            end
        end
        function set.align_signals(C, value)
            C.align_signals = value;
            C.D.auto_clip(true)
        end
    end
end
