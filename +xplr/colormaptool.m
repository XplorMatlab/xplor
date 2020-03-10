classdef colormaptool < xplr.graphnode
    
    properties (SetObservable=true)
        cmapdef
    end
    properties (SetObservable=true, AbortSet, SetAccess='private')
        cmap
    end
    properties (SetAccess='private')
        menu
    end
    
    events
        ChangedColorMap
    end
   
    methods
        function C = colormaptool(D)
            % Build menu
            buildmenu(C,D);
            % Set to default: jet
            C.cmap = 'jet';
        end
        function delete(C)
            delete@xplr.graphnode(C)
            if ~isprop(C,'menu'), return, end
            deleteValid(C.menu)
        end
        function buildmenu(C,D)
            hf = D.V.hf;
            
            % Create menu
            if isempty(C.menu)
                m = uimenu('parent',hf,'label','Color Map');
                C.menu = m;
            else
                m = C.menu;
                if nargin>=2 && get(m,'parent')~=hf, error 'existing menu does not have the same parent figure', end
                delete(get(m,'children'))
            end
            
            % Menu items
            mapnames = {'gray' 'jet' 'parula' 'mapgeog' 'mapgeogclip' 'mapclip' 'mapcliphigh' ...
                'mapcliplow' 'vdaq' 'green' 'red' 'blue-yellow' 'blue-red' 'maporient'};
            fn_propcontrol(C,'cmapdef',['menugroup' mapnames 'user...'],m);
            
            % Control visibility depending on dislay mode
            set(C.menu,'visible',fn_switch(D.displaymode,'image','on','off'));
            connectlistener(D,C,'displaymode','PostSet', ...
                @(u,e)set(C.menu,'enable',fn_switch(D.displaymode,'image','on','off')));
        end
    end
    
    % Color map and color map name: be careful with infinite loops!
    methods
        function set.cmapdef(C,x)
            % Handle colormap names and 'user..' command
            if ischar(x)
                name = x;
                switch name
                    case 'user...'
                        answer = inputdlg('Define color map here:','',1,{'rand(256,3)'});
                        if isempty(answer), return, end
                        try
                            x = evalin('base',answer{1});
                        catch
                            errordlg('Could not evaluate command to a valid color map','','modal')
                            return
                        end
                    otherwise
                        try
                            x = feval(name,256);
                        catch
                            errordlg(sprintf('''%s'' is not the name of a valid color map',name),'','modal')
                            return
                        end
                end
            else
                name = 'user';
            end
            % Set color map 
            if ~(isnumeric(x) && ismatrix(x) && size(x,2)==3 && all(x(:)>=0 & x(:)<=1))
                error('not a valid color map definition')
            end
            C.cmap = x; %#ok<MCSUP>
            % Set name
            C.cmapdef = name;
            % Notify
            notify(C,'ChangedColorMap')            
        end
    end
end