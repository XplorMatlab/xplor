classdef view < hgsetget
    % Main window
    %
    
    % Components 
    properties (SetAccess = 'private')
        hf              % figure
        slicer          % operation to get 'slice' from 'data'
        C               % control of data operation
        D               % main display
        panels          % panels for: display, control, lists, + some buttons
        context         % object that creates context menus on-the-fly
    end
    properties (SetObservable = true)
        controlvisible = false;  % logical - are the controls visible
    end
    
    % Dependent properties
    properties (Dependent, SetAccess = 'private')
        data
        slice
    end
        
    % Constructor
    methods
        function V = view(data,varargin)
            % DATA
            % define headers if needed
            if ~isa(data,'xplr.xdata')
                head = xplr.editHeader(data);
                if isempty(head)
                    % header editing was canceled
                    % (make Dependent properties 'data' and 'slice' valid
                    % to avoid errors when V will be displayed)
                    V.slicer = struct('data',[],'slice',[]);
                    return
                end
                data = xplr.xdata(data,head);
            end
            % create slicer
            if ~isa(data,'xplr.xdata'), error 'data argument must be a xplr.xdata object', end
            V.slicer = xplr.slicer(data);
            
            % LINKS
            % register view object to the bank
            xplr.bank.registerView(V)

            % PANELS
            % open figure and create panels
            init_panels(V)
            
            % CONTROL
            V.C = xplr.viewcontrol(V);
            
            % DISPLAY
            V.D = xplr.viewdisplay(V);
            
            % CONTEXT MENU FACILITY
            V.context = xplr.contextmenu(V);
            
%             % default filters
%             head = data.header;
%             F = xplr.xfilter.empty(1,0);
%             P = xplr.xpoint.empty(1,0);
%             for i=1:length(head)
%                 [Fi Pi] = xplr.bank.getFilterAndPoint(head(i));
% %                 % default selection: all
% %                 Fi.setSelection(xplr.selectionnd('all1D'))
%                 F = [F Fi]; %#ok<AGROW>
%                 P = [P Pi]; %#ok<AGROW>
%             end
%             V.slicer.setFilters(1:data.nd,F)
%             V.slicer.setPoints(1:data.nd,P)
%             
%             % settings relevant for the initialization?
%             opt = struct('in',[]);
%             Fopt = fieldnames(opt);
%             F = varargin(1:2:end);
%             idx = find(ismember(F,Fopt));
%             for i=idx
%                 f = F{i};
%                 opt.(f) = varargin{2*i};
%             end
%             varargin([idx idx+1]) = [];
%             
%             % init display
%             init_figure(V,opt.in)
%             init_axes(V)
%             init_pointer(V)
%             
%             % other settings
%             for k=1:2:length(varargin)
%                 set(V,varargin{k},varargin{k+1})
%             end

            % save object in base workspace
            assignin('base','V',V)
        end
        function delete(V)
            if ~isprop(V,'hf'), return, end
            deleteValid(V.hf,V.C,V.D)
        end
    end
    
    % Reset display
    methods
        function reset_display(V)
            delete(allchild(V.panels.display))
            V.D = xplr.viewdisplay(V);
        end
    end
    
    % Panels
    methods
        function init_panels(V)
            % figure
            V.hf = figure('integerhandle','off','handlevisibility','off','visible','off', ...
                'numbertitle','off','name','XPLOR','tag','XPLOR', ...
                'menubar','none', ...
                'color','white');
            delete(findall(V.hf,'parent',V.hf))
            set(V.hf,'DeleteFcn',@(u,e)delete(V))
            set(V.hf,'WindowButtonMotionFcn',@(u,e)donothing())
            if ~isempty(V.data.name), set(V.hf,'name',['XPLOR: ' V.data.name]), end
            
            % by defaults, callbacks should not be interruptible
            set(V.hf,'Interruptible','off', ...
                'DefaultUicontrolInterruptible','off', ...
                'DefaultLineInterruptible','off', ...
                'DefaultTextInterruptible','off')
            
            % panels: use panelorganizer tool
            V.panels.main = panelorganizer(V.hf,'H',2,[0 1],[0 1]); % left part has fixed width, right part adjusts automatically
            set(V.panels.main,'bordermode','push')
            V.panels.display = V.panels.main.setSubPanel(2);
            V.panels.allcontrols = V.panels.main.setSubOrg(1,'V',2,[1 1],[1 0]);
            connectlistener(V.panels.allcontrols.hobj,V,'Visible','PostSet',@(u,e)set(V,'controlvisible',get(V.panels.allcontrols.hobj,'visible')));
            V.panels.controlswidth = 120;
            V.panels.control = V.panels.allcontrols.setSubPanel(1);
            V.panels.listcombo = V.panels.allcontrols.setSubOrg(2,'H');
            V.controlvisible = true; % start with controls visible; it is necessary to explicitely set property to true (rather than having true as default value) to update display correctly
            
            % switch of control visibility
            ppi = get(0,'screenpixelsPerInch');
            controlon = fn_propcontrol(V,'controlvisible',{'pushbutton' {true false} {'>> Controls' '<< Controls'}}, ...
            	'parent',V.panels.display,'backgroundcolor','w','fontsize',800/ppi);
            fn_controlpositions(controlon.hu,V.panels.display,[0 1],[2 -17 65 15])
            
            % now figure can be visible
            set(V.hf,'visible','on')
        end
        function set.controlvisible(V,b)
            if strcmp(b,'toggle')
                b = ~V.controlvisible;
            else
                b = fn_switch(b,'logical');
                if b==V.controlvisible, return, end
            end
            % set property
            V.controlvisible = b;
            % positions
            main = V.panels.main;
            if b && main.extents(1)==0
                main.pushExtent(1,V.panels.controlswidth,'figleft'); %#ok<*MCSUP>
            elseif ~b && main.extents(1)>0
                V.panels.controlswidth = main.extents(1); % memorize for next controls re-opening
                main.pushExtent(1,0,'figleft');
            end
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
            if strcmp(get(in,'type'),'axes')
                V.ha = in;
                V.hf = get(V.ha,'parent');
            else
                fn_isfigurehandle(in) % create figure if in is a positive integer but not yet a figure handle
                V.hf = in;
                V.ha = axes('parent',V.hf);
            end
        end
        function init_axes(V)
            % structure with all secondary graphic objects
            s = struct;
            
            % axes on the side serve to generate "outside" mouse events
            col = get(V.hf,'color')*.9;
            s.ha2 = axes('parent',V.hf,'buttondownfcn',@(u,e)Mouse(V,'outsidedown'));
            fn_controlpositions(s.ha2(1),V.ha,[0 0 1 0], [0 -20 0 20]);
            s.ha2(2) = axes('parent',V.hf,'buttondownfcn',@(u,e)Mouse(V,'outsideleft'));
            fn_controlpositions(s.ha2(2),V.ha,[0 0 0 1], [-20 0 20 0]);
            s.ha2(3) = axes('parent',V.hf,'buttondownfcn',@(u,e)Mouse(V,'outsideboth'));
            fn_controlpositions(s.ha2(3),V.ha,[0 0 0 0], [-20 -20 20 20]);
            set(s.ha2,'handlevisibility','off','color',col, ...
                'xtick',[],'xcolor',col,'ytick',[],'ycolor',col)
            uistack(s.ha2,'bottom')
            
            % store structure
            V.hobj = s;
        end
        function init_pointer(V)
            h(1) = line('Parent',V.ha,'xdata',[0 0]);
            h(2) = line('Parent',V.ha,'ydata',[0 0]);
            h(3) = line('Parent',V.ha,'xdata',0,'ydata',0,'marker','.','linestyle','none'); % a single point
            V.hobj.pointer = h;
            if V.isimage
                set(h,'Color','white')
            else
                set(h,'Color','k')
                set(h(2:3),'visible','off')
            end
        end
    end
    
    
end

%---
function donothing()
% setting this function as main figure WindowButtonMotionFcn forces the
% figure to update CurrentPoint whenever the mouse is moved, and therefore
% have this property set currently even when clicking active controls
% (unfortunately, clicking a control does not set this property)
end
