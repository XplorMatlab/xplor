classdef headerEdit < hgsetget
    
    properties
        % data (will be read-only)
        dat
        sz
        nd
        
        % header under construction
        curhead = struct('sublabels',cell(1,0),'label',[], ...
            'unit',[],'start',[],'scale',[],'values',[],'colors',[], ...
            'isvalid',[],'markguess',[],'allguess',[]);
        header
        
        % graphics
        hf
        table
        ok
        uconfirm
        contextmenu
    end
    
    % Constructor
    methods
        function E = headerEdit(data)
            % data
            if nargin<1, data = rand(50,40,100); end
            if isa(data,'xplr.xdata')
                inputhead = data.header;
                E.dat = data.data;
                E.sz = data.sz;
            else
                inputhead = [];
                E.dat = data;
                E.sz = xplr.strictsize(data);
            end
            E.nd = length(E.sz);
            
            % guess header and store in a structure that admits invalid
            % definitions (units defined without start and scale being
            % defined)
            for i=1:E.nd
                % guesses or other start value
                if isempty(inputhead)
                    candidates = xplr.bank.getrecentheaders(E.sz(i));
                    okguess = ~isempty(candidates);
                    if ~okguess
                        candidates = xplr.header('',E.sz(i));
                    end
                else
                    okguess = false;
                    candidates = inputhead(i);
                end
                % store in a different format more adequate for display
                clear candidates2
                for j=1:length(candidates)
                    headj = candidates(j);
                    sublabels = {headj.sublabels.label};
                    values = headj.values; colors = [];
                    kcolor = strcmp(sublabels,'ViewColor');
                    if any(kcolor)
                        colors = cell2mat(values(:,kcolor));
                        sublabels(kcolor) = [];
                        values(:,kcolor) = [];
                    end
                    if isscalar(sublabels) && strcmp(sublabels{1},headj.label)
                        % show sublabels only if they are different from the
                        % summary label
                        sublabels = [];
                    end
                    candidates2(j) = struct('sublabels',{sublabels},'label',headj.label, ...
                        'unit',headj.unit,'start',headj.start,'scale',headj.scale, ...
                        'values',{values},'colors',colors, ...
                        'isvalid',true,'markguess',okguess,'allguess',[]); %#ok<AGROW>
                end
                % current head value
                occurence = sum(E.sz(1:i)==E.sz(i));
                if length(candidates2)>=occurence, j=occurence; else j=1; end
                E.curhead(i) = candidates2(j);
                if okguess
                    E.curhead(i).allguess = candidates2; % store all other alternatives as well
                end
            end
            
            E.hf = figure('integerhandle','off','handlevisibility','off', ...
                'numbertitle','off','name','Set headers', ...
                'menubar','none','resize','off');
            set(E.hf,'WindowButtonMotionFcn',@donothing) % force update of CurrentPoint when moving the mouse
            
            % menu for editing units
            menu_units(E)
            
            % init table
            init_table(E)
        end
    end
    
    % Edit measures
    methods
        function menu_units(E)
            delete(findall(E.hf,'type','uimenu','label','Edit recognized units'))
            m = uimenu('parent',E.hf,'label','Edit recognized units');
            measures = xplr.bank.getMeasures();
            labels = {measures.label};
            for i=1:length(labels)
                label = measures(i).label;
                units = {measures(i).units.unit};
                uimenu(m,'label',[label ' (' fn_strcat(units,',') ')'], ...
                    'callback',@(u,e)editMeasure(E,label))
            end
            uimenu(m,'label','Create a new set','separator',fn_switch(~isempty(labels)), ...
                'callback',@(u,e)editMeasure(E,''))
        end
        % this function is actually completely independent from E, except
        % it re-displays the menus at the end
        function editMeasure(E,oldlabel)
            
            % which measure to edit / new unit
            createnew = isempty(oldlabel);
            if createnew
                units = struct('unit',cell(1,0),'value',cell(1,0));
            else
                measures = xplr.bank.getMeasures();
                idx = strcmp(oldlabel,{measures.label});
                units = measures(idx).units;
            end
            
            % figure
            hfm = figure('windowstyle','modal');
            clf(hfm,'reset')
            set(hfm,'numbertitle','off','name','Edit Units', ...
                'menubar','none','resize','off')
            
            % some constants for control positioning
            dx = 10;
            dy = 10;
            h = 14;
            htab = 250;
            W = 220;
            H = dy+h+2*dy+htab+2*dy+2*h+dy;
            ycur = H;
            
            % use a fixed size for easy positioning of elements
            fn_setfigsize(hfm,[W H]);
            set(hfm,'defaultuicontrolunit','pixel')
            
            % label
            ycur = ycur-dy-h;
            uicontrol('style','text','string','Measure label', ...
                'pos',[dx ycur (W-2*dx)*.4 h]);
            xlabel = uicontrol('style','edit','string',oldlabel, ...
                'pos',[dx+(W-2*dx)*.4 ycur (W-2*dx)*.6 h]);
            if createnew, set(xlabel,'ToolTipString','e.g. time'), end
            
            % data (add an empty line)
            data = [column({units.unit}) column({units.value})];
            
            % table
            ycur = ycur-2*dy-htab;
            t = uitable('pos',[dx ycur W-2*dx htab]);
            set(t,'ColumnName',{'Unit' 'Value'},'ColumnFormat',{'char' 'numeric'},'ColumnEditable',true, ...
                'RowName','','ColumnWidth',repmat({(W-2*dx)/2-1},1,2),'Data',data)
            if createnew, set(t,'ToolTipString','e.g. unit: ms, value: 1e-3'), end
            set(t,'CellEditCallback',@checknewline), checknewline
            function checknewline(~,~)
                data = get(t,'Data');
                if isempty(data) || fn_find(data(end,:),'any')
                    data = [data; repmat({'' []},2,1)];
                    set(t,'Data',data)
                end
            end
            
            % ok button
            ycur = ycur-2*dy-2*h;
            okmeas = uicontrol('string','ok','callback',@(u,e)delete(u), ...
                'pos',[W-dx-50 ycur 50 2*h]);
            
            %             % make units of controls 'normalized' for meaningful figure
            %             % resize behavior
            %             set(get(hfm,'children'),'unit','normalized')
            
            % wait for use finished and grab data before closing window
            waitfor(okmeas)
            if ~ishandle(hfm), return, end % figure was closed, cancel
            label = get(xlabel,'String');
            data = get(t,'Data');
            close(hfm)
            
            % save edited measure in bank
            while ~isempty(data) && ~fn_find(data(end,:),'any'), data(end,:)=[]; end % remove empty lines
            if isempty(label) || isempty(data), return, end
            units = struct('unit',data(:,1),'value',data(:,2));
            if createnew
                xplr.bank.addMeasure(label,units);
            else
                xplr.bank.editMeasure(oldlabel,label,units);
            end
            
            % reinit menu to take into account the change in measures
            menu_units(E)
        end
    end
    
    % Edit headers
    methods
        function init_table(E)
            % figure
            W = 560; H = 380;
            fn_setfigsize(E.hf,W,H)
            
            % ok button
            E.ok = uicontrol('parent',E.hf,'string','ok','pos',[W-80 1 80 30], ...
                'callback',@(u,e)done(E));
            
            % "confirm all" button masks ok button! (indeed, ok will not be
            % enabled as long as there will remain dimensions to confirm)
            E.uconfirm = uicontrol('parent',E.hf,'string','Confirm all','callback',@(u,e)confirmall(E), ...
                'pos',[W-80 1 80 30]);
            uicontrol('parent',E.hf,'string','Reset all','callback',@(u,e)resetall(E), ...
                'pos',[W-2*80 1 80 30]);
            
            % empty table
            E.table = uitable('parent',E.hf,'pos',[1 31 W H-30], ...
                'ColumnName',{'Dim','Size','Label','Unit','Scale/Values','Colors',''}, ...
                'ColumnFormat',{'numeric' 'numeric' 'char' 'char' 'char' 'char' 'char'}, ...
                'ColumnEditable',logical([0 0 1 1 1 1 0]));
            set(E.table,'TooltipString','test')
            u = E.table;
            p = get(u,'pos'); w = p(3);
            widths = {30 55 70 70 [] 55 55};
            wavail = w - sum([widths{:}]) - 2;
            idxauto = find(fn_isemptyc(widths));
            [widths{idxauto}] = deal(floor(wavail/length(idxauto)));
            set(u,'ColumnWidth',widths)
            set(u,'RowName',[])
            set(u,'CellEditCallback',@(u,e)celledit(E,e), ...
                'CellSelectionCallback',@(u,e)cellselect(E,e))
            
            % fill table
            display_header(E,1:E.nd)
        end
        function display_header(E,idx)
            tdata = get(E.table,'Data');
            [iL iU iV iA iC] = columnIndices;
            for i = idx
                head = E.curhead(i);
                % dimension number
                tdata{i,1} = i;
                % length
                tdata{i,2} = E.sz(i);
                % header info
                [tdata{i,[iL iU iV iC]}] = display_headerinfo(head);
                % is guess?
                if length(head.allguess)>1
                    tdata{i,iA} = 'Choose...';
                elseif head.markguess 
                    tdata{i,iA} = 'Confirm';
                else
                    tdata{i,iA} = '';
                end            
            end
            set(E.table,'Data',tdata)
            % which values are guessed and need being confirmed
            color_table(E)
        end
        function color_table(E)
            % color rows with values that need being confirmed
            col = zeros(E.nd,3); allok = true;
            for i=1:E.nd
                if E.curhead(i).markguess
                    col(i,:) = [1 1 0]; allok = false;
                elseif ~E.curhead(i).isvalid
                    col(i,:) = [1 .4 .4]; allok = false;
                elseif mod(i,2)
                    col(i,:) = [1 1 1];
                else
                    col(i,:) = [1 1 1]*.94;
                end
            end
            if fn_matlabversion('newgraphics')
                set(E.table,'BackgroundColor',col)
            end
            
            % enabling of 'confirm all' and 'ok' buttons
            if ishandle(E.uconfirm) && ~any(any([E.curhead.markguess]))
                delete(E.uconfirm)
            end
            set(E.ok,'enable',fn_switch(allok))
        end
        function celledit(E,e)
            i = e.Indices(1);
            cnames = get(E.table,'ColumnName');
            ename = cnames{e.Indices(2)};

            % can remove the 'Confirm' flag
            E.curhead(i).markguess = false;
            
            % update header
            switch ename
                case 'Label'
                    update_label(E,i)
                case 'Unit'
                    update_unit(E,i)
                case 'Scale/Values'
                    update_value(E,i)
                case 'Colors'
                    update_colors(E,i)
            end
            
            % update display
            check_valid(E,i)
            display_header(E,i)
        end
        function update_label(E,i)
            head = E.curhead(i);
            tdata = get(E.table,'Data');
            iL = columnIndices;
            str = tdata{i,iL};
            tokens = regexp(str,'^(.*[^ ]) *\((.*)\)$','tokens');
            if ~isempty(tokens)
                % form 'label (label1*label2)'
                [head.label sublabels] = deal(tokens{1}{:});
                head.sublabels = fn_strcut(sublabels,'*');
            elseif any(str=='*')
                head.label = str;
                head.sublabels = fn_strcut(str,'*');
            else
                head.label = str;
                head.sublabels = [];
            end
            E.curhead(i) = head;
        end
        function update_unit(E,i)
            tdata = get(E.table,'Data');
            [~, iU] = columnIndices;
            unit = tdata{i,iU};
            tokens = regexp(unit,'^(.*) \[(.*)\]$','tokens');
            if ~isempty(tokens), unit = tokens{1}{1}; end
            head = E.curhead(i);
            % update start, scale, values if necessary
            if isempty(head.unit) && ~isempty(unit)
                [head.start head.scale head.values] = deal(1,1,{});
            elseif ~isempty(head.unit) && isempty(unit)
                [head.start head.scale head.values] = deal([],[],{});
            end
            head.unit = unit;
            E.curhead(i) = head;
        end
        function update_value(E,i)
            head = E.curhead(i);
            tdata = get(E.table,'Data');
            [~, ~, iV] = columnIndices;
            [type value] = read_value(tdata{i,iV},E.sz(i));
            switch type
                case 'invalid'
                    % do not accept change
                case 'measure'
                    [head.scale head.start] = dealc(value);
                    head.values = [];
                case 'enum'
                    [head.scale head.start head.values] = deal([]);
                case 'categorical'
                    [head.scale head.start] = deal([]);
                    head.values = value;
            end
            E.curhead(i) = head;
        end
        function update_colors(E,i)
            head = E.curhead(i);
            tdata = get(E.table,'Data');
            [~, ~, ~, ~, iC] = columnIndices;
            % colors not allowed for 'measure' header
            if ~isempty(head.unit), return, end
            % no colors?
            if isempty(tdata{i,iC}), E.curhead(i).colors = []; end
            % try setting colors
            try %#ok<TRYNC>
                colors = evalin('base',tdata{i,iC});
                if isequal(size(colors),[E.sz(i) 3]), E.curhead(i).colors = colors; end
            end
        end
        function check_valid(E,i)
            head = E.curhead(i);
            % number of labels
            [nlabel ncolumn] = deal(length(head.sublabels),size(head.values,2));
            if nlabel==0
                okv = (ncolumn<=1);
            else
                okv = (ncolumn==nlabel);
            end
            % if unit is defined, so must be start and scale
            okv = okv && (isempty(head.unit) || (~isempty(head.start) && ~isempty(head.scale)));
            E.curhead(i).isvalid = okv;
        end
        function cellselect(E,e)
            tdata = get(E.table,'Data');
            [iL iU iV iA iC] = columnIndices;
            if size(e.Indices,1)~=1, return, end
            i = e.Indices(1);
            headi = E.curhead(i);
            if ~isscalar(headi.allguess) && e.Indices(2)==iA
                % select among the list of all guesses: build and show a
                % menu with all possibilities
                deleteValid(E.contextmenu)
                m = uicontextmenu(E.hf);
                E.contextmenu = m;
                nguess = length(headi.allguess);
                for j=1:nguess
                    [label unit scale_value color] = display_headerinfo(headi.allguess(j));
                    lab = fn_strcat({label unit scale_value color},'; ');
                    uimenu(m,'label',lab,'callback',@(u,e)useguess(E,i,j))
                end
                set(m,'pos',get(E.hf,'CurrentPoint'),'visible','on')
            elseif headi.markguess
                if e.Indices(2)==iA
                    % confirm
                    tdata{i,iA} = '';
                    set(E.table,'Data',tdata)
                    color_table(E)
                elseif any(e.Indices(2)==[iL iU iV iC])
                    % reset
                    [tdata{i,[iL iU iV iC iA]}] = deal('');
                    [headi.label headi.unit headi.start headi.scale headi.values ...
                        headi.colors] = deal('','',[],[],{},[]);
                end
                set(E.table,'Data',tdata)
                headi.markguess = false;
                headi.allguess = [];
                E.curhead(i) = headi;
                color_table(E)
            end
        end
        function useguess(E,i,j)
            % update header
            allguess = E.curhead(i).allguess;
            E.curhead(i) = allguess(j);
            E.curhead(i).markguess = false; % do not mark as guess, since the user specifically selected it
            E.curhead(i).allguess = []; % it would also make sense not to remove 'allguess' and thus leave the possibility to the user to change again to another guess
            % update display
            display_header(E,i)
        end
        function confirmall(E)
            tdata = get(E.table,'Data');
            [~, ~, ~, iA] = columnIndices;
            tdata(:,iA) = {''};
            set(E.table,'Data',tdata)
            [E.curhead.markguess] = deal(false);
            color_table(E)
        end
        function resetall(E)
            tdata = get(E.table,'Data');
            [iL iU iV iA] = columnIndices;
            tdata(:,[iL iU iV iA]) = {''};
            set(E.table,'Data',tdata)
            [E.curhead.markguess] = deal(false);
            [E.curhead.label] = deal('');
            [E.curhead.unit] = deal('');
            [E.curhead.start] = deal([]);
            [E.curhead.scale] = deal([]);
            for i=1:E.nd
                E.curhead(i).values = cell(E.sz(i),0);
            end
            color_table(E)
        end
        function done(E)
            % build headers
            E.header = xplr.header.empty(1,0);
            for i=1:E.nd
                head = E.curhead(i);
                if isempty(head.unit)
                    % categorical
                    if isempty(head.sublabels), head.sublabels={head.label}; end
                    if ~isempty(head.colors)
                        head.sublabels{end+1} = 'ViewColor';
                        head.values(:,end+1) = num2cell(head.colors,2);
                    end
                    if isempty(head.values)
                        E.header(i) = xplr.header(head.label,E.sz(i));
                    else
                        E.header(i) = xplr.header(head.label,head.sublabels,head.values);
                    end
                else
                    % measure
                    [unit, ~, measure] = read_unit(head.unit);
                    if isempty(measure)
                        dimlabel = xplr.dimensionlabel(head.label,'numeric',unit);
                    else
                        dimlabel = xplr.dimensionlabel(head.label,'numeric',measure.units);
                    end
                    E.header(i) = xplr.header(dimlabel,E.sz(i),head.start,head.scale);
                end
            end
            
            assignin('base','ans',E.header)
            evalin('base','display(ans)')
            
            % close figure -> calling editHeader function can proceed
            delete(E.hf)
            drawnow
        end
    end
    
end

%---
function donothing(u,e)
% this function is set to the WindowButtonMotionFcn property of the figure
% to force update of CurrentPoint when moving the cursor
end

%---
function [iL iU iV iA iC] = columnIndices

[iL iU iV iA iC] = deal(3,4,5,7,6);

end

%---
function [unit str measure] = read_unit(unit)
% if unit is recognized as a unit for a known measure, return an enhanced
% string display and the list of all units for this measure

tokens = regexp(unit,'^(.*) \[(.*)\]$','tokens');
if ~isempty(tokens), unit = tokens{1}{1}; end

[measurelabel, ~, measure] = xplr.bank.getunitinfo(unit);
isknown = ~isempty(measurelabel);
if isknown
    comment = [' [' measurelabel ']'];
else
    comment = [];
end
str = [unit comment];

end

%---
function [type value] = read_value(scale_value,n)

% empty?
if isempty(scale_value)
    type = 'enum';
    value = [];
    return
end

% scale + start?
tokens = regexp(scale_value,'^(.*)\[ *(start){0,1}(.*)\]$','tokens');
if ~isempty(tokens)
    tok = tokens{1};
    try %#ok<TRYNC>
        scale = evalin('base',tok{1});
        start = evalin('base',tok{3});
        type = 'measure';
        value = [scale start];
        return
    end
end

% try evaluating in base workspace
try x = evalin('base',scale_value); catch, x = []; end % can also be a string that evaluates to a number

% number?
if isscalar(x)
    type = 'measure';
    value = [x 0];
    return
elseif isvector(x) && isnumeric(x) && length(x)==n && max(abs(diff(x,2)))<diff(x(1:2))/1e6
    % vector of values
    if max(abs(diff(x,2)))==0
        scale = diff(x(1:2));
    else
        scale = (x(n)-x(1))/(n-1);
    end
    start = x(1) - scale;
    type = 'measure';
    value = [scale start];
    return
end

% list of items?
if isvector(x) && n>1, x = column(x); end
if iscell(x) && size(x,1)==n
    list = x;
else
    sep = fn_switch(any(scale_value==','),',',' ');
    list = fn_strcut(scale_value,sep);
    if isvector(list) && n>1, list = column(list); end
end
if length(list)==n
    type = 'categorical';
    value = list;
    return
end

% failed to read string
type = 'invalid';
value = [];
return

end

%---
function [label unit scale_value color] = display_headerinfo(head)

% label
if isempty(head.sublabels)
    label = head.label;
else
    autolabel = fn_strcat(head.sublabels,'*');
    if strcmp(autolabel,head.label)
        label = head.label;
    else
        label = [head.label ' (' autolabel ')'];
    end
end
% unit
[~, unit] = read_unit(head.unit);
% scale/values
if isempty(head.scale)
    scale_value = display_value('categorical',head.values);
else
    scale_value = display_value('measure',[head.scale head.start]);
end
% color
if isempty(head.colors)
    color = '';
else
    color = sprintf('[%ix%i array]',size(head.colors));
end

end

%---
function str = display_value(type,value)

switch type
    case 'measure'
        str = [num2str(value(1),12) ' [start ' num2str(value(2),12) ']'];
    case 'categorical'
        if isempty(value)
            str = '';
        elseif isvector(value)
            str = fn_strcat(value,',');
        else
            str = sprintf('%ix%i table',size(value));
        end
end

end