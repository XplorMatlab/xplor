classdef contextmenu < hgsetget
% A context menu that is created only when it is being raised
    
    properties (Access='private')
        V       % parent view
        menu    % context menu handle
    end
    
    % Constructor and on-the-fly creation
    methods
        function M = contextmenu(V)
            % contextmenu constructor
            % contextmenu (V), take the view as argument set it as parent
            % then set matlab uicontextmenu 
            
            M.V = V;
            M.menu = uicontextmenu('parent',V.hf);
        end
        function raise(M,flag,varargin)
        % raise(M,flag,varargin) displays the uicontextmenu
        % M is the contextmenu
        % flag can be 'label' if it's the context menu of the label of a dimension of the graph
        % or 'datadim' if it's the context menu of the dimensions list next to the graph 
        % varargin
            
            % Hide and delete previous entries of this context menu
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
            ZS = D.zoomslicer;
            
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
                    
                    % select ZoomFilter key (check the created menu item
                    % that corresponds to the current key)
                    m2 = uimenu(m,'label','zoom filter','Separator',fn_switch(docolor));
                    availablekeys = xplr.bank.availableFilterKeys();
                    newkey = max(availablekeys)+1;
                    keyvalues = [0 availablekeys newkey];
                    fn_num2str(availablekeys, 'shared zoom %i', 'cell');
                    keydisplays = [ ...
                        'private zoom' ...
                        fn_num2str(availablekeys, 'shared zoom %i', 'cell') ...
                        num2str(newkey,'shared zoom %i (new key)')
                        ];
                    curkey = D.zoomfilters(dim).linkkey;
                    for i=1:length(keyvalues)
                        keyvalue = keyvalues(i);      
                        uimenu(m2,'label',keydisplays{i},'checked',fn_switch(isequal(curkey,keyvalue)), ...
                            'callback',@(u,e)ZS.changeKey(dim,keyvalue));
                    end
                    
                    % select crossSelector key 
                    curfilt = D.navigation.dimfilters{dim};
                    if ~isempty(curfilt)
                        m2 = uimenu(m,'label','cross selector key','Separator',fn_switch(docolor));
                        
                        availablekeys = xplr.bank.availableFilterKeys();
                        newkey = max(availablekeys)+1;
                        keyvalues = [0 availablekeys newkey];
                        fn_num2str(availablekeys, 'cross selector key %i', 'cell');
                        keydisplays = [ ...
                            'private cross selector' ...
                            fn_num2str(availablekeys, 'cross selector key %i', 'cell') ...
                            num2str(newkey,'cross selector key %i (new key)')
                            ];
                            curkey = curfilt.linkkey;
                        uimenu(m2, 'label', 'show point selector','callback',@(u,e)dimaction(D.V.C,'showFilterPointWindow',curkey,dim));
                        for i=1:length(keyvalues)
                            keyvalue = keyvalues(i);
                            uimenu(m2,'label',keydisplays{i},'checked',fn_switch(isequal(curkey,keyvalue)), ...
                                'callback',@(u,e)connectPointFilter(D.navigation,dim,keyvalue));
                        end
                    end
                        
                case 'datadim'
                    % Change filters
                    
                    
                    % display the available keys to apply a new or existing
                    % filter. 
                    m2 = uimenu(m,'label','Add/Change filter');
                    uimenu(m,'label','Remove Filters','callback',@(u,e)dimaction(C,'rmfilter',1,dim))
                    availablekeys = xplr.bank.availableFilterKeys();
                    newkey = max(availablekeys)+1;
                    keyvalues = [0 availablekeys newkey];
                    fn_num2str(availablekeys, 'shared zoom %i', 'cell');
                    keydisplays = [ ...
                        'private filter' ...
                        fn_num2str(availablekeys, 'shared filter %i', 'cell') ...
                        num2str(newkey,'shared filter %i (new key)')
                        ];
                    for i=1:length(keyvalues)
                        keyvalue = keyvalues(i);      
                        uimenu(m2,'label',keydisplays{i}, ...
                            'callback',@(u,e)dimaction(C,'filter',keyvalue,dim));
                    end
                    
                    % uimenu(m,'label','Synchronize Zoom','callback',)
                    % m1 = uimenu(m,'label','scroll wheel','separator','on');
                    %    L.menuitems.scrollwheel = uimenu(m1,'label','activated', ...
                    %'checked',L.scrollwheel,'callback',@(u,e)set(L,'scrollwheel',fn_switch(L.scrollwheel,'toggle')));
                    %    uimenu(m1,'label','make default in figure', ...
                    % 'callback',@(u,e)set(L,'scrollwheel','default'));
                    
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
