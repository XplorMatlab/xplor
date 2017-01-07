classdef cliptool < hgsetget
    
    % properties controlled by (and as seen from) the menu
    properties (SetObservable=true, AbortSet, SetAccess='private')
        autoclipmode = 'minmax';
    end
    properties (SetObservable=true, AbortSet)
        autoclipmodeNocenter = 'minmax'
        center = []             % [], 0 or 1
        adjust = 'mean(line)'   % 'none', 'mean' or 'mean(line)'
        span = 'local'          % 'curview', 'local' or 'shared (sharename)'
    end
    % other properties
    properties (SetAccess='private')
        sharename = ''
        menu
    end
    
    events
        ChangedClip
    end
   
    methods
        function C = cliptool(hf)
            if nargin==0
                hf = figure(826);
                clf(hf)
                set(hf,'handlevisibility','off')
            end
            buildmenu(C,hf)
        end
        function delete(C)
            if ~isprop(C,'menu'), return, end
            deleteValid(C.menu)
        end
        function buildmenu(C,hf)
            if isempty(C.menu)
                m = uimenu('parent',hf,'label','Clipping');
                C.menu = m;
            else
                m = C.menu;
                if nargin>=2 && get(m,'parent')~=hf, error 'existing menu does not have the same parent figure', end
                delete(get(m,'children'))
            end
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
            % set clip
            uimenu(m,'label','Set Clip Manually...','callback',@(u,e)setManualclip(C))
            uimenu(m,'label','Do Auto-Clip','callback',@(u,e)notify(C,'ChangedClip',xplr.eventinfo('clip','automode')))
            % adjust
            fn_propcontrol(C,'adjust', ...
                {'menuval' {'none' 'mean' 'mean(line)'} {'None' 'Adjust each Element by its Mean Value' 'Auto (Adjust Time Courses but not Images)'} {'None' 'Mean' 'Auto'}}, ...
                m,'label','Adjust','separator','on');
            % span
            fn_propcontrol(C,'span', ...
                {'menuval' {'curview' 'local' 'chgshare'} {'Local, Adjust to Current View' 'Local' 'Shared...'} {'Current View' 'Local' ''}}, ...
                m,'label','Span');
        end
    end
    
    % Setting clip manually
    methods
        function setManualclip(C,clip)
            if nargin<2
                clip = fn_input('Enter min and max',[0 1]);
                if length(clip)~=2, return, end % value not valid, do not use it
            end
            notify(C,'ChangedClip',xplr.eventinfo('clip','clip',clip))
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
            % notify the change
            notify(C,'ChangedClip',xplr.eventinfo('clip','automode'))
        end
        function set.center(C,center)
            % set property
            C.center = center;
            % update autoclip mode
            C.autoclipmode = C.autoclipmodeNocenter;
            if ~isempty(C.center), C.autoclipmode = [C.autoclipmode '[' num2str(C.center) ']']; end
            % notify the change
            notify(C,'ChangedClip',xplr.eventinfo('clip','automode'))
        end
    end
    
    % Adjusting
    methods
        function set.adjust(C,adj)
            C.adjust = adj;
            notify(C,'ChangedClip',xplr.eventinfo('clip','adjust'))
        end
    end
    
    % Sharing
    methods
        function set.span(C,span)
            if strcmp(span,'chgshare')
                C.sharename = fn_input('Clip sharing name','#1');
                C.span = ['Shared (' C.sharename ')']; %#ok<*MCSUP>
            else
                C.span = span;
            end
            notify(C,'ChangedClip',xplr.eventinfo('clip','span'))
            % shared clipping not really handled yet
        end
    end
end