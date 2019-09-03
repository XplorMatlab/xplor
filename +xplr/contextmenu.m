classdef contextmenu < hgsetget
    % A context menu that is created only when it is being raised
    
    properties (Access='private')
        % parent view
        V
        % context menu handle
        menu
    end
    
    % Constructor and on-the-fly creation
    methods
        function M = contextmenu(V)
            M.V = V;
            M.menu = uicontextmenu('parent',V.hf);
        end
        function raise(M,flag,varargin)
            % Delete previous entries
            set(M.menu,'visible','off')
            delete(get(M.menu,'children'))
            % Create the entries
            createentries(M,flag,varargin{:})
            % Show the menu
            p = get(M.V.hf,'CurrentPoint');
            set(M.menu,'pos',p,'visible','on')
        end
    end
    methods (Access='private')
        function createentries(M,flag,dim)
            % Shortcuts
            m = M.menu;
            C = M.V.C;
            D = M.V.D;
            
            % Dimension
            head = D.zslice.header(dim);
            
            % Create entries
            switch flag
                case 'label'
                    % Line properties
                    % (color)
                    docolor = strcmp(D.displaymode,'time courses');
                    if docolor
                        uimenu(m,'label',['Color according to ' head.label],'checked',fn_switch(isequal(D.colordim,dim)), ...
                            'callback',@(u,e)D.setColorDim(fn_switch(isequal(D.colordim,dim),[],dim)))
                        uimenu(m,'label','Display color legend', ...
                            'enable',fn_switch(isequal(D.colordim,dim)),'checked',fn_switch(D.showcolorlegend), ...
                            'callback',@(u,e)set(D,'showcolorlegend',~D.showcolorlegend))
                    end
                    
                    % Binning
                    m1 = uimenu(m,'label','Binning','Separator',fn_switch(docolor));
                    binvalues = {1 2 3 4 'set'};
                    bindisplays = {'none' '2' '3' '4' 'other...'};
                    curbin = D.zoomfilters(dim).bin;
                    for i=1:length(binvalues)
                        bin = binvalues{i};
                        uimenu(m1,'label',bindisplays{i},'checked',fn_switch(isequal(curbin,bin)), ...
                            'callback',@(u,e)setbin(D,dim,bin));
                    end
                case 'datadim'
                    % Change filters
                    uimenu(m,'label','Add/Show Filters','callback',@(u,e)dimaction(C,'filter',1,dim))
                    uimenu(m,'label','Remove Filters','callback',@(u,e)dimaction(C,'rmfilter',1,dim))
                    uimenu(m,'label','Add private filter','callback',@(u,e)dimaction(C,'filter',0,dim))
                    % uimenu(m,'label','Synchronize Zoom','callback',)
                    %m1 = uimenu(m,'label','scroll wheel','separator','on');
                    %    L.menuitems.scrollwheel = uimenu(m1,'label','activated', ...
                    %'checked',L.scrollwheel,'callback',@(u,e)set(L,'scrollwheel',fn_switch(L.scrollwheel,'toggle')));
                    %    uimenu(m1,'label','make default in figure', ...
                    %'callback',@(u,e)set(L,'scrollwheel','default'));
                otherwise
                    error('unknown flag ''%s''',flag)
            end
            
        end
    end
    
end


%---
function setbin(D,d,bin)

if strcmp(bin,'set')
    bin = fn_input('Binning',D.zoomfilters(d).bin,'stepper 1 1 Inf 1');
    if isempty(bin), return, end
end
D.zoomfilters(d).setBin(bin)

end


