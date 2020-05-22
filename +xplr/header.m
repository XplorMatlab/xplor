classdef header < hgsetget
    % function H = header(label[,unit],n,start,scale)   [measure]
    % function H = header(label,unit,n[,start,scale])   [measure]
    % function H = header(label,n)                      [categorical]
    % function H = header([label,]sublabels,values)     [categorical]
    % function H = header({args header 1},{args header 2},...)
    %---
    % Container for header information in a single dimension.
    %
    % Input:
    % * label   a dimensionlabel object or a string (in the latter case
    %           a dimensionlabel object will be created)
    % * unit    unit; add a '+' sign (e.g. 'ms+') to try recognizing the
    %           unit and storing other available units for this measure
    % * n       expected number of data elements
    % * start   numerical value for first data element
    % * scale   numerical value for step between data elements
    % * values  a cell array with n rows and as many columns as there are
    %           labels
    %
    % There are two types of header; this type is automatically inferred
    % from the syntax of the call to 'header' constructor.
    % * 'measure' header typically corresponds to a continuous dimension
    %   along which data is acquired at regularly-spaced intervals (e.g.
    %   time or space), though it can be discrete as well (e.g. day number)
    % * 'categorical' header corresponds to samples without regular
    %   organizations; information about this samples can be given using
    %   a table of values that can be multi-column (e.g. date and location)
    % 
    % A header object, despite being a handle, cannot be modified by the
    % user (this would be too risky, as the same header can be shared by
    % multiple data or filter objects). Instead of modifying a header,
    % it should be replaced by a newly created header.
    
    properties (SetAccess='private')
        sublabels = xplr.dimensionlabel.empty(1,0);      % a dimensionlabel object
        label           % char
        n               % length
        categorical     % boolean
        start           % scalar - for measure only
        scale           % scalar - for measure only
        values          % cell array (categorical) or vector (measure)
    end
    properties (Dependent, Transient)
        unit
        allunits
        type
        ismeasure
        isenum
        iscategoricalwithvalues
        ncolumn
    end
    % The properties below are computed on the fly and then stored
    properties (Access='private')
        ID              % access with method getID
        measureSpaceID  % access with metho getMeasureSpaceID
        itemnames       % access with method getItemNames
    end
    
    % Constructor, copy, display
    methods
        function H = header(varargin)
            % special cases
            if nargin==0        
                % empty input -> empty header (nd=1 !)
                return
            elseif nargin==1 && isnumeric(varargin{1})
                n = varargin{1};
                if n>1
                    H(n) = xplr.header();
                end
                return
            elseif nargin==1 && isa(varargin{1},'xplr.header')
                H1 = varargin{1};
                H = xplr.header(length(H1));
                H.copyin(H1)
                return
            end
            
            % multiple header definitions?
            if ~iscell(varargin{1})
                ismultiple = false;
            elseif nargin~=2
                ismultiple = true;
            else
                % we must still distinguish between the syntax for 2
                % headers, or a single categorical header
                ismultiple = ~all(fn_map(@ischar,varargin{1}));
            end
            if ismultiple
                for i=1:nargin
                    H(i) = xplr.header(varargin{i}{:}); %#ok<AGROW>
                end
                return
            end
                
            % categorical?
            H.categorical = (nargin<3) || (nargin==3 && ~isnumeric(varargin{3}));
            if H.categorical
                % categorical
                % input: is a global name given? is table empty?
                if nargin==2
                    [lab, table_or_n] = deal(varargin{:});
                    if ischar(lab) && isscalar(table_or_n) && isnumeric(table_or_n)
                        % only a number of samples, no table description
                        table = cell(table_or_n,0);
                        name = lab;
                        lab = {};
                    else
                        table = table_or_n;
                        if isnumeric(table), table = num2cell(table); end
                        name = [];
                    end
                elseif nargin==3
                    [name, lab, table] = deal(varargin{:});
                else
                    error 'not enough input arguments'
                end
                % check size of table and assign to 'values' property
                if ischar(lab), lab = {lab}; end
                nlabel = length(lab);
                if (nlabel==1 && ~isvector(table)) || (nlabel>1 && size(table,2)~=nlabel)
                    error 'values should be a cell array with as many columns as there are labels'
                elseif nlabel==1
                    table = column(table);
                end
%                 if nlabel==0 % if table is empty, create one column with label 'index'
%                     lab = {'index'};
%                     table = num2cell((1:size(table,1))');
%                 end
                H.values = table;
                H.n = size(table,1);
                % set sublabels
                if isa(lab,'xplr.dimensionlabel')
                    H.sublabels = lab;
                else
                    % need to infer the class of each label
                    if H.n==0 && nlabel>0
                        error 'no elements, cannot infer the type of label(s), please specify directly by using dimensionlabel object(s)'
                    end
                    for i=1:nlabel
                        xi = table{1,i};
                        type = xplr.dimensionlabel.infertype(xi);
                        H.sublabels(i) = xplr.dimensionlabel(lab{i},type);
                    end
                end
                % set summary label
                if ~isempty(name)
                    H.label = name;
                else
                    H.label = fn_strcat({H.sublabels.label},'*');
                end
            else
                % measure
                % set unique sublabel and label
                lab = varargin{1};
                if isa(lab,'xplr.dimensionlabel')
                    if ~strcmp(lab.type,'numeric')
                        error 'type of dimensionlabel must be ''numeric'' for non-categorical header'
                    end
                    H.sublabels = lab;
                else
                    % is a unit given?
                    if ischar(varargin{2})
                        unit = varargin{2};
                        varargin(2) = [];
                        if length(unit)>1 && unit(end)=='+'
                            % check the bank for other linked units
                            unit(end) = [];
                            [~, ~, measure] = xplr.bank.getunitinfo(unit);
                            if ~isempty(measure), unit = measure.units; end
                        end
                        H.sublabels = xplr.dimensionlabel(lab,'numeric',unit);
                    else
                        H.sublabels = xplr.dimensionlabel(lab,'numeric');
                    end
                end
                H.label = H.sublabels.label;
                % set number of samples, start and scale
                switch length(varargin)
                    case 2
                        H.n = varargin{2};
                        [H.start H.scale] = deal(1);
                    case 4
                        [H.n H.start H.scale] = deal(varargin{2:4});
                    otherwise
                        error 'incorrect header specification'
                end
%                 % set values
%                 H.values = num2cell(H.start + (0:H.n-1)'*H.scale);
            end
        end
        function disp(H)
            valid = ~isempty(H(1).n);
            if ~valid
                disp@hgsetget(H)
                return
            end
            str = [class(H) ' object with following info:'];
            if ~isscalar(H)
                sz = size(H);
                str = [fn_strcat(sz,'x') ' ' str];
            end
            fprintf(['  ' str '\n\n'])
            fields = {'label' 'type' 'n' 'unit'};
            if isscalar(H)
                if H.ismeasure
                    fields = {'label' 'type' 'n' 'unit' 'allunits' 'start' 'scale'};
                else
                    fields = {'label' 'type' 'n' 'sublabels' 'values'};
                end
            end
            nf = length(fields);
            val = cell(1,nf);
            for i=1:nf
                f = fields{i};
                switch f
                    case 'allunits'
                        if H.ismeasure
                            if isempty(H.allunits)
                                val{i} = '<invalid empty allunits>'; % old version of xplr.header class
                            else
                                val{i} = fn_strcat({H.allunits.unit},','); 
                            end
                        end
                        if isempty(val{i}), val{i} = '[]'; end
                    case 'sublabels'
                        val{i} = fn_strcat({H.sublabels.label},',');
                    case 'values'
                        okval = false;
                        if H.ncolumn==1 && H.n<20
                            try val{i} = fn_strcat(H.values,','); okval = true; end %#ok<TRYNC>
                        end
                        if ~okval
                            val{i} = sprintf('%ix%i cell array',size(H.values));
                        end
                    otherwise
                        str = fn_strcat({H.(f)},',');
                        if isempty(str), str = '[]'; end
                        val{i} = str;
                end
            end
            names = fliplr(char(fn_map(@fliplr,fields)));
            for i=1:nf
                disp(['    ' names(i,:) ': ' val{i}])
            end
            fprintf('\n')
        end
    end
    methods (Access='protected')
        function H1 = copy(H)
            H1 = xplr.header;
            H1.copyin(H);
        end
        function copyin(H1,H)
            [H1.sublabels] = deal(H.sublabels);
            [H1.label] = deal(H.label);
            [H1.n] = deal(H.n);
            [H1.categorical] = deal(H.categorical);
            [H1.start] = deal(H.start);
            [H1.scale] = deal(H.scale);
            [H1.values] = deal(H.values);
            % the copy is intended to be followed by some modifications,
            % therefore, do not copy ID and itemnames, which will become
            % invalid when these modifications will occur
        end
    end
    
    % Dependent and computed properties
    methods
        function u = get.unit(H)
            if H.categorical
                u = '';
            else
                u = H.sublabels.unit;
            end
        end
        function u = get.allunits(H)
            if H.ismeasure
                u = H.sublabels.allunits;
            else
                u = [];
            end
        end
        function type = get.type(H)
            if H.categorical
                type = 'categorical';
            else
                type = 'measure';
            end
        end
        function b = get.ismeasure(H)
            b = ~H.categorical;
        end
        function b = get.isenum(H)
            b = H.categorical && (H.ncolumn==0);
        end
        function b = get.iscategoricalwithvalues(H)
            b = H.categorical && (H.ncolumn>0);
        end
        function n = get.ncolumn(H)
            n = size(H.values,2);
        end
    end
    
    % Header comparisons
    methods
        function b = isequal(H1,H2)
            % cannot use the default isqual function, because property
            % itemnames or ID might be computed for one header and not for
            % the other
            if ~isscalar(H1) || ~isscalar(H2)
                b = false;
                if ~isequal(size(H1),size(H2)), return, end
                for i=1:numel(b), if ~isequal(H1(i),H2(i)), return, end, end
                b = true;
                return
            end
            b = (H1.n==H2.n) && isequal(H1.values,H2.values) ... % start with the most likely to be unequal
                && isequal(H1.label,H2.label) && isequal(H1.sublabels,H2.sublabels) ...
                && isequal({H1.categorical,H1.start,H1.scale},{H2.categorical,H2.start,H2.scale});
        end
        function b = is_equal(H1,H2)
            % function b = is_equal(H1,H2)
            %---
            % alias to isequal, needed when both H1 and H2 are instances of
            % subclasses of xplr.header
            b = isequal(H1,H2);
        end
        function ID = getID(H)
            % Get a unique identifier that identifies the header (we will
            % have H1 == H2 if and only if H1.getID() == H2.getID())
            if ~isscalar(H)
                ID = zeros(size(H));
                for i=1:numel(H), ID(i) = getID(H(i)); end
                return
            end
            if isempty(H.ID)
                H.ID = fn_hash({H.n,H.sublabels,H.label,H.categorical,H.start,H.scale,H.values},'num'); 
            end
            ID = H.ID;
        end
        function ID = getMeasureSpaceID(H)
            % Get a unique identifier that identifies the space inside
            % which the data dimensions described by the header lies: this
            % is a hash number of the header's label and unit.
            if ~isscalar(H)
                ID = zeros(size(H));
                for i=1:numel(H), ID(i) = getMeasureSpaceID(H(i)); end
                return
            end
            if H.ismeasure
                if isempty(H.measureSpaceID)
                    H.measureSpaceID = fn_hash({H.label,H.unit},'num'); 
                end
                ID = H.measureSpaceID;
            else
                ID = [];
            end
        end
    end
    
    % Group measure headers operating on the same space
    methods
        function connections = measure_grouping(H)
            % function connections = measure_grouping(H)
            %---
            % returns square boolean matrix indicating with ones the
            % dimensions whose headers are measure with the same units
            nh = length(H);
            connections = false(nh,nh);
            idxmeasure = find([H.ismeasure]);
            for i = idxmeasure
                for j = setdiff(idxmeasure,i)
                    connections(i,j) = strcmp(H(i).unit, H(j).unit);
                end
            end
        end
    end
    
    % Access value in table
    methods
        function x = getValue(H,label,idx)
            % function x = getValue(H,label[,idx])
            % access value(s) in the table (categorical header only)
            if H.ismeasure, error 'measure header does not have sample values', end
            icol = strcmp({H.sublabels.label},label);
            if ~any(icol)
                disp(sprintf('header has no sub-label ''%s''',label)) %#ok<DSPS>
                x = [];
                return
            end
            if nargin<3
                x = H.values(:,icol);
            else
                if ~isscalar(idx), error 'index must be scalar', end
                x = H.values{idx,icol};
            end
        end
    end
    
    % List of labels for each item
    methods
        function str = getItemNames(H,idx)
            doall = (nargin<2);
            if doall, idx = 1:H.n; end
            if isempty(H.itemnames)
                % compute names
                nval = length(idx);
                if H.ncolumn>0
                    % any column 'name' in the values table?
                    idxname = find(strcmpi({H.sublabels.label},'name'));
                    if isempty(idxname), idxname = 1; end
                    itemvalues = H.values(:,idxname);
                    str = cell(1,nval);
                    for k=1:nval
                        val = itemvalues{idx(k)};
                        if isempty(val)
                            str{k} = num2str(idx(k));
                        elseif ischar(val)
                            str{k} = val;
                        elseif isnumeric(val) || islogical(val)
                            str{k} = fn_idx2str(val,':,');
                            if length(str{k})>12, str{k} = [str{k}(1:10) '...']; end
                        elseif iscell(val)
                            str{k} = fn_strcat(val,',');
                        else
                            error 'cannot form string from value'
                        end
                    end
                elseif H.categorical
                    str = fn_num2str(idx,'cell')';
                else % (measure)
                    str = fn_num2str(H.start+(idx-1)*H.scale,['%.4g' H.unit],'cell')';
                end
                if doall, H.itemnames = str; end
            elseif doall
                str = H.itemnames;
            else
                str = H.itemnames(idx);
            end
        end
    end
    
    % Header update
    methods
        function newhead = updateHeader(H,flag,ind,value)
            % function newhead = updateHeader(H,flag,ind,value)
            if ~fn_ismemberstr(flag,{'all' 'new' 'chg' 'remove' 'chg&new' 'chg&rm' 'perm' 'chgdim'})
                error('cannot update header for flag ''%s''',flag)
            end
            newhead = copy(H);
            switch H.type
                case 'measure'
                    error 'not implemented yet (call updateMeasureUpdate?)'
                case 'categorical'
                    switch flag
                        case {'all' 'chgdim'}
                            newhead.n = length(ind);
                            if H.ncolumn==0
                                newhead.values = cell(newhead.n,0);
                            else
                                if isvector(value) && length(value)==newhead.n
                                    value = column(value);
                                elseif size(value,1)~=newhead.n
                                    error 'number of rows of header values does not match header length'
                                elseif size(value,2)~=newhead.ncolumn
                                    error 'number of columns of header values does not match number of labels'
                                end
                                newhead.values = value;
                            end
                        case 'new'
                            newhead.n = H.n + length(ind);
                            if H.ncolumn==0
                                newhead.values = cell(newhead.n,0);
                            else
                                newhead.values(ind,:) = value;
                            end
                        case 'chg'
                            % no change in size, change values only if
                            % specified
                            if nargin>=4, newhead.values(ind,:) = value; end
                        case 'chg&new'
                            newhead.n = H.n + length(ind{2});
                            if H.ncolumn==0
                                newhead.values = cell(newhead.n,0);
                            else
                                newhead.values([ind{:}],:) = value;
                            end
                        case 'chg&rm'
                            newhead.n = H.n - length(ind{2});
                            if nargin>=4, newhead.values(ind{1},:) = value; end
                            newhead.values(ind{2},:) = [];
                        case 'remove'
                            newhead.n = H.n - length(ind);
                            newhead.values(ind,:) = [];
                        case 'perm'
                            newhead.values = H.values(ind,:);
                    end
            end
        end
        function newhead = updateMeasureHeader(H,newn,newstart,newscale)
            % function newhead = updateMeasureHeader(H,newn[,newstart[,newscale]])
            newhead = copy(H);
            newhead.n = newn;
            if nargin>=3, newhead.start = newstart; end
            if nargin>=4, newhead.scale = newscale; end
        end
        function checkHeaderUpdate(H,flag,ind,newhead)
            % check that the new header seems to meet with the specified
            % change
            
            % flag 'chgdim' indicates a complete change, i.e. there is just
            % nothing to be checked
            if strcmp(flag,'chgdim'), return, end
            
            % first check that sublabels and type are the same (newhead
            % however might have more or less sublabels)
            if xor(newhead.categorical,H.categorical)
                error 'new header is not of the same type as current header'
            end
            nsub = min(length(H.sublabels), length(newhead.sublabels));
            if ~isequal(newhead.sublabels(1:nsub),H.sublabels(1:nsub))
                error 'new header has different sublabel(s) from current header'
            end
            
            % that's it for 'all' flag
            if strcmp(flag,'all'), return, end
            
            % check number of points
            switch flag
                case 'new'
                    ncheck = H.n + length(ind);
                case 'remove'
                    ncheck = H.n - length(ind);
                case 'chg&new'
                    ncheck = H.n + length(ind{2});
                case 'chg&rm'
                    ncheck = H.n - length(ind{2});
                otherwise
                    ncheck = H.n;
            end
            if newhead.n~=ncheck
                error 'wrong number of elements in new header'
            end
            
            % check values
            if H.ismeasure
                switch flag
                    case 'all'
                        % nothing more needs to be checked
                    case 'chg'
                        % change is possible only if all values were
                        % changed
                        ok = (newhead.start==H.start && newhead.scale==H.scale) || (length(ind)==H.n);
                    case {'chg&new' 'new' 'chgdata'}
                        ok = (newhead.start==H.start) && (newhead.scale==H.scale);
                    case {'remove' 'chg&rm'}
                        if strcmp(flag,'remove'), idxrm = ind; else idxrm = ind{2}; end
                        ok = (newhead.scale==H.scale);
                        if ok
                            if idxrm(1)==1 && idxrm(end)==H.n
                                nrmstart = find(diff(idxrm)~=1);
                                if isempty(nrmstart)
                                    % removed all!
                                    ok = (newhead.start==H.start);
                                elseif isscalar(nrmstart)
                                    ok = abs(newhead.start-(H.start+nrmstart*H.scale))<abs(H.scale)/1000;
                                else
                                    ok = false;
                                end
                            elseif idxrm(1)==1
                                ok = abs(newhead.start-(H.start+length(idxrm)*H.scale))<abs(H.scale)/1000 && all(diff(idxrm)==1);
                            elseif idxrm(end)==H.n
                                ok = (newhead.start==H.start) && all(diff(idxrm)==1);
                            else
                                ok = false;
                            end
                        end
                    case 'perm'
                        ok = false;
                    otherwise
                        error('invalid flag ''%s'' for header update',flag)
                end
            else
                switch flag
                    case 'all'
                        % nothing more needs to be checked
                        [idxH idxN] = deal([]);
                    case 'chg'
                        unchg=true(1,H.n); unchg(ind)=false;
                        [idxH idxN] = deal(find(unchg)); % (much) faster than setdiff
                    case {'new' 'chgdata'}
                        [idxH idxN] = deal(1:H.n);
                    case 'chg&new'
                        unchg=true(1,H.n); unchg(ind{1})=false;
                        [idxH idxN] = deal(find(unchg));
                    case 'chg&rm'
                        unchg=true(1,H.n); unchg([ind{:}])=false;
                        idxH = find(unchg);
                        unchg=true(1,newhead.n); unchg([ind{1}])=false;
                        idxN = find(unchg);
                    case 'remove'
                        unchg=true(1,H.n); unchg(ind)=false;
                        idxH = find(unchg);
                        idxN = 1:newhead.n;
                    case 'perm'
                        idxH = ind;
                        idxN = 1:H.n;
                    otherwise
                        error('invalid flag ''%s'' for header update',flag)
                end
                ok = isequal(newhead.values(idxN,1:H.ncolumn),H.values(idxH,:));
            end
            if ~ok, error 'New header values are not consistent with the operation', end
        end
        function newhead = addLabel(H,label,value)
            % xplr.dimensionlabel object
            if isa(label,'xplr.dimensionlabel')
                name = label.label;
                if nargin<3, value = label.defaultval; end
            else
                name = label;
                if iscell(value), x = value{1}; else x = value; end % there must be a value to infer type correctly
                type = xplr.dimensionlabel.infertype(x);
                label = xplr.dimensionlabel(name,type);
            end
            if any(strcmp(name,{H.sublabels.label}))
                error('a sub-label ''%s'' already exists',name)
            end
            % add it to the list of labels
            idx = length(H.sublabels)+1;
            newhead = copy(H);
            newhead.sublabels(idx) = label;
            if iscell(value)
                newhead.values(:,idx) = value;
            else
                if size(newhead.values,1)==0
                    newhead.values = cell(0,idx);
                else
                    [newhead.values{:,idx}] = deal(value);
                end
            end
        end
        function newvalues = trackValues(H,indices)
            nind = length(indices);
            newvalues = cell(nind,H.ncolumn);
            for i=1:nind
                idx = indices{i};
                v = H.values(idx,:);
                if isscalar(idx)
                    % just copy values
                    newvalues(i,1:H.ncolumn) = v;
                else
                    % several indices together
                    for j=1:H.ncolumn
                        switch H.sublabels(j).label
                            case 'ViewColor'
                                % average colors
                                newvalues{i,j} = mean(cat(1,v{:,j}),1);
                            case 'Name'
                            otherwise
                                % default behavior: compute union of values
                                switch H.sublabels(j).type
                                    case {'numeric' 'logical'}
                                        newvalues{i,j} = unique([v{:,j}],'stable');
                                    case 'char'
                                        newvalues{i,j} = unique(v(:,j)','stable');
                                end
                        end
                    end
                end
            end
        end
    end
    
end