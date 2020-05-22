classdef cliptool < xplr.graphnode
% cliptool 
    
    % properties controlled by (and as seen from) the menu
    properties (SetAccess='private')
        D
    end
    properties (SetObservable=true, AbortSet, SetAccess='private')
        autoclipmode = 'minmax'
    end
    properties (SetObservable=true, AbortSet)
        autoclipmodeNocenter = 'minmax'     % 
        center = []                         % [], 0 or 1
        independent_dimID_mem = []          % dimension ID of dimensions along which clipping is not uniform
        align_signals = ''                  % '', 'nmean', 'nmedian'
        adjust_to_view = false
    end
    properties (Dependent, SetAccess='private')
        independent_dim
        linked_dim
    end

    properties (SetAccess='private')
        sharename = ''                      % 
        menu                                % 
    end
    
    % Constructor, menu
    methods
        function C = cliptool(D)
            % constructor
            C.D = D;
            C.menu = uimenu('parent',D.V.hf,'label','Clipping', ...
                'callback',@(u,e)clip_menu(C));
        end
        function delete(C)
            delete@xplr.graphnode(C)
            if ~isprop(C,'menu'), return, end
            deleteValid(C.menu)
        end
        function clip_menu(C)
            m = C.menu;
            delete(get(m,'children'))
            % auto-clip specifications
            % (create a submenu whose label will update automatically)
            P = fn_propcontrol(C,'autoclipmode', ...
                {'menuval' {}}, ...
                m,'label','Auto-Clip Method');
            m1 = P.hparent;
            % (some possible values)
            fn_propcontrol(C,'autoclipmodeNocenter', ...
                {'menugroup' {'minmax' 'std1' 'std2' 'std3' 'std5' 'std.5' 'std.2' 'prc.1' 'prc1' 'prc5'} ...
                {'Min to Max' '1 STD' '2 STD' '3 STD' '5 STD' '1/2 STD' '1/5 STD' '[.1% 99.9%]' '[1% 99%]' '[5% 95%]'}}, ...
                m1);
            % (all possibilities through additional input dialogs)
            uimenu(m1,'label','Use Mean and STD...','separator','on','callback',@(u,e)setAutoclipMode(C,'setstd'))
            uimenu(m1,'label','Use Percentiles...','callback',@(u,e)setAutoclipMode(C,'setprc'))
            % (fix the center)
            fn_propcontrol(C,'center', ...
                {'menugroup' {0 1 []} {'center on 0' 'center on 1'}}, ...
                m1);
            % auto-clip
            if ~isempty(C.independent_dim)
                uimenu(m,'label','Do Auto-Clip (current cell(s) only)','callback',@(u,e)C.D.autoClip(false))
                uimenu(m,'label','Do Auto-Clip (all cells)','callback',@(u,e)C.D.autoClip(true))
            else
                uimenu(m,'label','Do Auto-Clip','callback',@(u,e)C.D.autoClip(true))
            end
            % adjust options
            fn_propcontrol(C,'adjust_to_view','menu', ...
                {'parent',m,'label','Always adjust clipping to current view','separator','on'});
            external_dim = C.D.external_dim;
            if ~isempty(external_dim)
                if isscalar(external_dim)
                    checked = ismember(external_dim,C.independent_dim);
                    uimenu(m,'label','Independent clipping range for each grid cell', ...
                        'checked',checked, ...
                        'callback',@(u,e)setIndependentDim(C,external_dim,~checked))
                else
                    m1 = uimenu('parent',m,'label','Independent clipping range for dimension(s)');
                    for dimID = external_dim
                        checked = ismember(dimID,C.independent_dim);
                        uimenu(m1,'label',C.D.slice.dimensionLabel(dimID), ...
                            'checked',checked, ...
                            'callback',@(u,e)setIndependentDim(C,dimID,~checked))
                    end
                    uimenu(m1,'label','(all)','separator','on', ...
                        'callback',@(u,e)setIndependentDim(C,external_dim,true))
                    uimenu(m1,'label','(none)', ...
                        'callback',@(u,e)setIndependentDim(C,external_dim,false))
                end
            end
            if strcmp(C.D.displaymode,'time courses')
                fn_propcontrol(C,'align_signals', ...
                    {'menuval' {'' 'nmean' 'nmedian'} {'(do not center)' 'on their mean' 'on their median'} {'(none)' 'mean' 'median'}}, ...
                    {'parent',m,'label','Align signals'});
            end
            % manual clip
            if ~isempty(C.independent_dim)
                uimenu(m,'label','Set Clip Manually... (current cell(s) only)','separator','on', ...
                    'callback',@(u,e)C.setManualclip([],false))
                uimenu(m,'label','Set Clip Manually... (all cells)', ...
                    'callback',@(u,e)C.setManualclip([],true))
            else
                uimenu(m,'label','Set Clip Manually...','separator','on', ...
                    'callback',@(u,e)C.setManualclip([],true))
            end
        end
    end
    
    % Get dependent
    methods
        function dim = get.independent_dim(C)
            indp_dim = C.D.slice.dimensionNumber(C.independent_dimID_mem);
            dim = intersect(C.D.external_dim, indp_dim);
        end
        function dim = get.linked_dim(C)
            indp_dim = C.D.slice.dimensionNumber(C.independent_dimID_mem);
            dim = setdiff(C.D.external_dim, indp_dim);
        end
    end
    
    % Setting clip manually
    methods
        function setManualclip(C,clip,all_cells)
            if nargin<2 || isempty(clip)
                clip = fn_input('Enter min and max',[0 1]);
                if length(clip)~=2, return, end % value not valid, do not use it
            end
            % update display
            if nargin<3, all_cells = true; end
            C.D.setClip(clip,all_cells)
        end
    end
    
    % Setting autoclip mode
    methods
        function setAutoclipMode(C,comp)
            switch comp
                case 'setstd'
                    nstd = fn_input('Mean +/- N*std',1,'stepper 1 0 Inf');
                    if isempty(nstd) || nstd==0, return, end
                    C.autoclipmodeNocenter = ['std' num2str(nstd)];
                case 'setprc'
                    prc = fn_input('Percentiles (one or two values between 0 and 100)',1);
                    if isempty(prc), return, end
                    switch length(prc)
                        case 1
                            C.autoclipmodeNocenter = ['prc' num2str(prc)];
                        case 2
                            C.autoclipmodeNocenter = ['prc' num2str(prc(1)) '-' num2str(prc(2))];
                        otherwise
                            % value not valid, do not set property
                    end
                otherwise
                    C.autoclipmodeNocenter = comp;
                    return
            end
        end
        function set.autoclipmodeNocenter(C,comp)
            % set property
            C.autoclipmodeNocenter = comp;
            % add the centering information to final autoclip mode
            C.autoclipmode = C.autoclipmodeNocenter;
            if ~isempty(C.center), C.autoclipmode = [C.autoclipmode '[' num2str(C.center) ']']; end
            % update display
            C.D.autoClip(true)
        end
        function set.center(C,center)
            % set property
            C.center = center;
            % update autoclip mode
            C.autoclipmode = C.autoclipmodeNocenter;
            if ~isempty(C.center), C.autoclipmode = [C.autoclipmode '[' num2str(C.center) ']']; end
            % update display
            C.D.autoClip(true)
        end
    end
    
    % Adjusting
    methods
        function setIndependentDim(C,dimID,independent)
            dimID = C.D.slice.dimensionID(dimID);
            if independent
                C.independent_dimID_mem = union(C.independent_dimID_mem,dimID);
            else
                C.independent_dimID_mem = setdiff(C.independent_dimID_mem,dimID);
            end
            % update display
            C.D.autoClip(true)
        end
        function set.adjust_to_view(C,value)
            C.adjust_to_view = value;
            % update display
            if value
                C.D.autoClip(true)
            end
        end
        function set.align_signals(C,value)
            C.align_signals = value;
            C.D.autoClip(true)
        end
    end
end