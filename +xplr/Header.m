classdef Header < handle
    % function H = header(label, n, unit, scale[, start])   [measure]
    % function H = header(label, unit, n, scale[, start])   [measure]
    % function H = header(dim_label, n, scale[, start])     [measure]
    % function H = header(label, n)                         [categorical]
    % function H = header([label,] sub_labels, values)      [categorical]
    % function H = header({args header 1}, {args header 2}, ...)
    % function H = header([n1 n2 ...], {simplified args header 1}, {simplified args header 2}, ...)
    %---
    % Container for header information in a single dimension.
    %
    % Input:
    % * label   a dimension_label object or a string (in the latter case
    %           a dimension_label object will be created)
    % * unit    unit; add a '+' sign (e.g. 'ms+') to try recognizing the
    %           unit and storing other available units for this measure
    % * n       expected number of data elements
    % * start   numerical value for first data element
    % * scale   numerical value for step between data elements
    % * values  a cell array with n rows and as many columns as there are
    %           labels
    % 
    % "simplified args header" refers to arguments as above but without 'n'
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
        sub_labels = xplr.DimensionLabel.empty(1, 0);      % a dimension_label object
        label           % char
        n               % length
        categorical     % boolean
        start           % scalar - for measure only
        scale           % scalar - for measure only
        values          % cell array (categorical) or vector (measure)
    end
    properties (Dependent, Transient)
        unit
        all_units
        type
        is_measure
        is_datetime
        is_enum
        is_categorical_with_values
        n_column
    end
    % The properties below are computed on the fly and then stored
    properties (Access='private')
        id              % access with method get_id
        measure_space_id  % access with metho get_measure_space_id
        item_names       % access with method get_item_names
    end
    
    % Constructor, copy, display
    methods
        function H = Header(varargin)
            % function H = header(label,n,unit,scale[,start])   [measure]
            % function H = header(label,unit,n,scale[,start])   [measure]
            % function H = header(dimlabel,n,scale[,start])     [measure]
            % function H = header(label,n)                      [categorical]
            % function H = header([label,]sublabels,values)     [categorical]
            % function H = header({args header 1},{args header 2},...)
            % function H = header()
            % function H = header(n)
            % function H2 = header(H1)

            % special cases
            if nargin == 0
                % empty input -> empty, invalid, header (but not empty
                % array of headers: use xplr.header.empty(size) to define
                % an empty array of headers)
                return
            elseif nargin == 1 && isnumeric(varargin{1})
                % array of empty headers
                n = varargin{1};
                if n>1
                    H(n) = xplr.Header();
                end
                return
            elseif nargin==1 && isa(varargin{1}, 'xplr.Header')
                % copy
                H1 = varargin{1};
                H = xplr.Header(length(H1));
                H.copy_in(H1)
                return
            end
            
            % multiple header definitions?
            % (special case: starting with size specifications)
            if isnumeric(varargin{1})
                % create xplr.header objects from labels
                spec = varargin(2:end);
                nd = length(spec);
                sz = varargin{1};
                assert(length(sz) <= nd, 'less header specifications than number of dimensions')
                sz(end+1:nd) = 1;
                for i=1:nd
                    if iscell(spec{i})
                        args = spec{i};
                        if ~ischar(args{end}) && ~isscalar(args{end})
                            % categorical header defined by a table of
                            % values: no need to add number of samples
                            H(i) = xplr.Header(args{:});
                        else
                            % other header formats: add number of samples
                            H(i) = xplr.Header(args{1}, sz(i), args{2:end});
                        end
                    else
                        label = spec{i};
                        assert(ischar(label), 'error argument')
                        H(i) = xplr.Header(label, sz(i));
                    end
                end
                return
            end
            % (regular case)
            if ~iscell(varargin{1})
                is_multiple = false;
            elseif nargin ~= 2
                is_multiple = true;
            else
                % we must still distinguish between the syntax for 2
                % headers, or a single categorical header
                is_multiple = ~all(brick.map(@ischar, varargin{1}));
            end
            if is_multiple
                for i = 1:nargin
                    H(i) = xplr.Header(varargin{i}{:}); %#ok<AGROW>
                end
                return
            end
                
            % categorical?
            if nargin < 3
                H.categorical = true;
            elseif nargin == 3
                num2 = isnumeric(varargin{2}) || isdatetime(varargin{2}) || isduration(varargin{2});
                num3 = isnumeric(varargin{2}) || isdatetime(varargin{2}) || isduration(varargin{2});
                H.categorical = (~num2 && ~num3);
            else
                H.categorical = false;
            end
            if H.categorical
                build_categorical_header(H,varargin{:})
            else
                build_measure_header(H,varargin{:})
            end
        end
        function disp(H)
            valid = ~isempty(H(1).n);
            if ~valid
                disp@handle(H)
                return
            end
            str = [class(H), ' object with following info:'];
            if ~isscalar(H)
                sz = size(H);
                str = [brick.strcat(sz,'x'), ' ', str];
            end
            fprintf(['  ', str, '\n\n'])
            fields = {'label', 'type', 'n', 'unit'};
            if isscalar(H)
                if H.is_measure
                    fields = {'label', 'type', 'n', 'unit', 'all_units', 'start', 'scale'};
                else
                    fields = {'label', 'type', 'n', 'sub_labels', 'values'};
                end
            end
            nf = length(fields);
            val = cell(1, nf);
            for i=1:nf
                f = fields{i};
                switch f
                    case 'all_units'
                        if H.is_measure
                            if isempty(H.all_units)
                                val{i} = '';
                            else
                                val{i} = brick.strcat({H.all_units.unit}, ',');
                            end
                        end
                        if isempty(val{i}), val{i} = '[]'; end
                    case 'sub_labels'
                        val{i} = brick.strcat({H.sub_labels.label}, ',');
                    case 'values'
                        ok_val = false;
                        if H.n_column==1 && H.n<20
                            try val{i} = brick.strcat(H.values,',');
                                ok_val = true;
                            end %#ok<TRYNC>
                        end
                        if ~ok_val
                            val{i} = sprintf('%ix%i cell array', size(H.values));
                        end
                    otherwise
                        str = brick.strcat({H.(f)}, ',');
                        if isempty(str), str = '[]'; end
                        val{i} = str;
                end
            end
            names = fliplr(char(brick.map(@fliplr, fields)));
            for i=1:nf
                disp(['    ', names(i,:), ': ', val{i}])
            end
            fprintf('\n')
        end
    end
    methods (Access='protected')
        function build_categorical_header(H, varargin)
            % function H = header(label, n)                      [categorical]
            % function H = header([label,] sublabels, values)    [categorical]

            % input: is a global name given? is table empty?
            if length(varargin)==2
                [lab, table_or_n] = deal(varargin{:});
                if ischar(lab) && isscalar(table_or_n) && isnumeric(table_or_n)
                    % only a number of samples, no table description
                    table = cell(table_or_n, 0);
                    name = lab;
                    lab = {};
                else
                    table = table_or_n;
                    if isnumeric(table) || isdatetime(table)
                        table = num2cell(table); 
                    elseif isstring(table)
                        table = cellstr(table);
                    elseif istable(table)
                        table = table2cell(table);
                    elseif ~iscell(table)
                        error('incorrect type ''%s'' for table', class(table))
                    end
                    name = [];
                end
            elseif length(varargin)==3
                [name, lab, table] = deal(varargin{:});
            else
                error 'not enough input arguments'
            end
            % check size of table and assign to 'values' property
            if ischar(lab), lab = {lab}; end
            n_label = length(lab);
            if (n_label==1 && ~isvector(table)) || (n_label>1 && size(table,2)~=n_label)
                error 'values should be a cell array with as many columns as there are labels'
            elseif n_label==1
                table = brick.column(table);
            end
            %                 if n_label==0 % if table is empty, create one column with label 'index'
            %                     lab = {'index'};
            %                     table = num2cell((1:size(table,1))');
            %                 end
            H.values = table;
            H.n = size(table,1);
            % set sublabels
            if isa(lab,'xplr.DimensionLabel')
                H.sub_labels = lab;
            else
                % need to infer the class of each label
                if H.n==0 && n_label>0
                    error 'no elements, cannot infer the type of label(s), please specify directly by using dimensionlabel object(s)'
                end
                for i=1:n_label
                    xi = table{1,i};
                    type = xplr.DimensionLabel.infer_type(xi);
                    H.sub_labels(i) = xplr.DimensionLabel(lab{i}, type);
                end
            end
            % set summary label
            if ~isempty(name)
                H.label = name;
            else
                H.label = brick.strcat({H.sub_labels.label}, '*');
            end
        end
        function build_measure_header(H,varargin)
            % function H = header(label, n, unit, scale[, start])   [measure]
            % function H = header(label, unit, n, scale[, start])   [measure]
            % function H = header(dimlabel, n, scale[, start])      [measure]

            % input -> set number of samples, scale and start
            lab = varargin{1};
            unit_ = [];
            H.scale = []; H.start = [];
            for i = 2:length(varargin)
                a = varargin{i};
                if ischar(a)
                    unit_ = a;
                elseif isempty(H.n)
                    H.n = a;
                elseif isempty(H.scale)
                    H.scale = a;
                else
                    H.start = a;
                end
            end
            if isempty(H.n), error 'number of samples not defined', end
            if isempty(H.scale), H.scale = 1; end
            if isempty(H.start), H.start = 0; end

            % set label, which is also the unique sublabel
            if isa(lab, 'xplr.DimensionLabel')
                assert(strcmp(lab.type,'numeric'), 'type of DimensionLabel must be ''numeric'' for non-categorical header')
                assert(isempty(unit_), 'if first input is a DimensionLabel object, unit must not be re-defined')
                H.sub_labels = lab;
            else
                % is a unit given?
                if ~isempty(unit_)
                    if length(unit_)>1 && unit_(end)=='+'
                        % check the bank for other linked units
                        unit_(end) = [];
                        [~, ~, measure] = xplr.Bank.get_unit_info(unit_);
                        if ~isempty(measure), unit_ = measure.units; end
                    end
                    H.sub_labels = xplr.DimensionLabel(lab, 'numeric', unit_);
                else
                    H.sub_labels = xplr.DimensionLabel(lab, 'numeric');
                end
            end
            H.label = H.sub_labels.label;
        end
        function H1 = copy(H)
            H1 = xplr.Header;
            H1.copy_in(H);
        end
        function copy_in(H1, H)
            [H1.sub_labels] = deal(H.sub_labels);
            [H1.label] = deal(H.label);
            [H1.n] = deal(H.n);
            [H1.categorical] = deal(H.categorical);
            [H1.start] = deal(H.start);
            [H1.scale] = deal(H.scale);
            [H1.values] = deal(H.values);
            % the copy is intended to be followed by some modifications,
            % therefore, do not copy id and item_names, which will become
            % invalid when these modifications will occur
        end
    end
    
    % Dependent and computed properties
    methods
        function u = get.unit(H)
            if H.categorical || H.is_datetime
                u = '';
            else
                u = H.sub_labels.unit;
            end
        end
        function u = get.all_units(H)
            if H.is_measure
                u = H.sub_labels.all_units;
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
        function b = get.is_measure(H)
            b = ~H.categorical;
        end
        function b = get.is_datetime(H)
            b = H.is_measure() && isdatetime(H.start);
        end
        function b = get.is_enum(H)
            b = H.categorical && (H.n_column == 0);
        end
        function b = get.is_categorical_with_values(H)
            b = H.categorical && (H.n_column > 0);
        end
        function n = get.n_column(H)
            n = size(H.values, 2);
        end
    end
    
    % Header comparisons
    methods
        function b = isequal(H1, H2)
            % cannot use the default isqual function, because property
            % item_names or id might be computed for one header and not for
            % the other
            if ~isscalar(H1) || ~isscalar(H2)
                b = false;
                if ~isequal(size(H1),size(H2)), return, end
                for i=1:numel(b), if ~isequal(H1(i), H2(i)), return, end, end
                b = true;
                return
            end
            b = (H1.n == H2.n) && isequal(H1.values, H2.values) ... % start with the most likely to be unequal
                && isequal(H1.label, H2.label) && isequal(H1.sub_labels, H2.sub_labels) ...
                && isequal({H1.categorical, H1.start, H1.scale}, {H2.categorical, H2.start, H2.scale});
        end
%        function b = isequal(H1, H2)
%            % function b = isequal(H1,H2)
%            %---
%            % alias to isequal, needed when both H1 and H2 are instances of
%            % subclasses of xplr.header
%            b = isequal(H1, H2);
%        end
        function id = get_id(H)
            % Get a unique identifier that identifies the header (we will
            % have H1 == H2 if and only if H1.get_id() == H2.get_id())
            if ~isscalar(H)
                id = zeros(size(H));
                for i=1:numel(H), id(i) = get_id(H(i)); end
                return
            end
            if isempty(H.id)
                H.id = brick.hash({H.n, H.sub_labels, H.label, H.categorical, H.start, H.scale, H.values}, 'num');
            end
            id = H.id;
        end
        function id = get_measure_space_id(H)
            % Get a unique identifier that identifies the space inside
            % which the data dimensions described by the header lies: this
            % is a hash number of the header's label and unit.
            if ~isscalar(H)
                if ~all([H.is_measure])
                    id = [];
                    return
                end
                id = zeros(size(H));
                for i=1:numel(H)
                    id(i) = get_measure_space_id(H(i));
                end
                return
            end
            if H.is_measure && ~H.is_datetime
                if isempty(H.measure_space_id)
                    H.measure_space_id = brick.hash({H.label,H.unit}, 'num');
                end
                id = H.measure_space_id;
            else
                id = [];
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
            idx_measure = find([H.is_measure] & [H.n]>1);
            for i = idx_measure
                for j = setdiff(idx_measure, i)
                    connections(i, j) = strcmp(H(i).unit, H(j).unit);
                end
            end
        end
    end
    
    % Access value in table
    methods
        function x = get_value(H, label, idx)
            % function x = get_value(H,label[,idx])
            % access value(s) in the table (categorical header only)
            if H.is_measure, error 'measure header does not have sample values', end
            icol = strcmp({H.sub_labels.label}, label);
            if ~any(icol)
                disp(sprintf('header has no sub-label ''%s''', label)) %#ok<DSPS>
                x = [];
                return
            end
            if nargin<3
                x = H.values(:, icol);
            else
                if ~isscalar(idx), error 'index must be scalar', end
                x = H.values{idx, icol};
            end
        end
        function c_map = get_color(H, idx)
            % function c_map = get_color(H[,idx])
            k_color = strcmp({H.sub_labels.label}, 'ViewColor');
            if nargin < 2
                idx = 1:H.n;
            end
            if any(k_color)
                c_map = cell2mat(H.values(idx, k_color));
            else
                c_map = brick.colorset('plot12', idx);
            end
        end
    end
    
    % List of labels for each item
    methods
        function names = get_item_names(H, idx)
            do_all = (nargin < 2);
            if do_all, idx = 1:H.n; end
            if isempty(H.item_names)
                % compute names
                n_val = length(idx);
                if H.n_column > 0
                    % any column 'name' in the values table?
                    idx_name = find(strcmpi({H.sub_labels.label}, 'name'));
                    if isempty(idx_name), idx_name = 1; end
                    item_values = H.values(:, idx_name);
                    if idx_name > 1
                        empty = brick.isemptyc(item_values);
                        item_values(empty) = H.values(empty, 1);
                    end
                    names = cell(1, n_val);
                    for k=1:n_val
                        val = item_values{idx(k)};
                        if isempty(val)
                            names{k} = num2str(idx(k));
                        elseif ischar(val)
                            names{k} = val;
                        elseif isnumeric(val) || islogical(val)
                            names{k} = brick.idx2str(val, ':,');
                            if length(names{k}) > 12, names{k} = [names{k}(1:10), '...']; end
                        elseif iscell(val)
                            names{k} = brick.strcat(val, ',');
                        elseif isdatetime(val)
                            names{k} = char(val);
                        else
                            error 'cannot form string from value'
                        end
                    end
                elseif H.categorical
                    names = brick.num2str(idx, 'cell')';
                elseif H.is_datetime
                    % (datetime -> guess what is the best format based on start and end)
                    v = H.start + (0:H.n-1) * H.scale;
                    step = H.scale;
                    v.Format = xplr.auto_datetime_format(v(1), v(end), step);
                    H.item_names = cellstr(char(v));
                    names = H.item_names(idx);
                else
                    % (measure)
                    names = brick.num2str(H.start + (idx-1)*H.scale, ['%.4g', H.unit], 'cell')';
                end
                if do_all, H.item_names = names; end
            elseif do_all
                names = H.item_names;
            else
                names = H.item_names(idx);
            end
        end
        function name = get_item_name(H,idx)
            if ~isscalar(idx)
                error 'call get_item_names method for multiple indices'
            end
            names = H.get_item_names(idx);
            name = names{1};
        end
    end
    
    % Header update
    methods
        function new_head = update_header(H, flag, ind, value)
            % function new_head = update_header(H,flag,ind,value)
            if strcmp(flag, 'sub_data')
                % header did not change
                new_head = H;
                return
            end
            if ~brick.ismemberstr(flag, {'all', 'new', 'chg', 'remove', 'chg&new', 'chg&rm', 'perm', 'chg_dim'})
                error('cannot update header for flag ''%s''', flag)
            end
            new_head = copy(H);
            switch H.type
                case 'measure'
                    % we do nothing, simply assuming that scale and start
                    % did not change!
                case 'categorical'
                    switch flag
                        case {'all', 'chg_dim'}
                            new_head.n = length(ind);
                            if H.n_column == 0
                                new_head.values = cell(new_head.n, 0);
                            else
                                if isvector(value) && length(value) == new_head.n
                                    value = brick.column(value);
                                elseif size(value, 1) ~= new_head.n
                                    error 'number of rows of header values does not match header length'
                                elseif size(value, 2) ~= new_head.n_column
                                    error 'number of columns of header values does not match number of labels'
                                end
                                new_head.values = value;
                            end
                        case 'new'
                            new_head.n = H.n + length(ind);
                            if H.n_column == 0
                                new_head.values = cell(new_head.n, 0);
                            else
                                new_head.values(ind, :) = value;
                            end
                        case 'chg'
                            % no change in size, change values only if
                            % specified
                            if nargin >= 4, new_head.values(ind, :) = value; end
                        case 'chg&new'
                            new_head.n = H.n + length(ind{2});
                            if H.n_column == 0
                                new_head.values = cell(new_head.n, 0);
                            else
                                new_head.values([ind{:}], :) = value;
                            end
                        case 'chg&rm'
                            new_head.n = H.n - length(ind{2});
                            if nargin >= 4, new_head.values(ind{1}, :) = value; end
                            new_head.values(ind{2}, :) = [];
                        case 'remove'
                            new_head.n = H.n - length(ind);
                            new_head.values(ind, :) = [];
                        case 'perm'
                            new_head.values = H.values(ind, :);
                    end
            end
        end
        function new_head = update_measure_header(H, new_n, new_start, new_scale)
            % function new_head = update_measure_header(H,new_n[,new_start[,new_scale]])
            new_head = copy(H);
            new_head.n = new_n;
            if nargin >= 3, new_head.start = new_start; end
            if nargin >= 4, new_head.scale = new_scale; end
        end
        function check_header_update(H, flag, ind, new_head)
            % check that the new header seems to meet with the specified
            % change
            
            % flag 'chg_dim' indicates a complete change, i.e. there is just
            % nothing to be checked
            if strcmp(flag, 'chg_dim'), return, end
            
            % first check that sub_labels and type are the same (new_head
            % however might have more or less sub_labels)
            if xor(new_head.categorical, H.categorical)
                error 'new header is not of the same type as current header'
            end
            nsub = min(length(H.sub_labels), length(new_head.sub_labels));
            if ~isequal(new_head.sub_labels(1:nsub), H.sub_labels(1:nsub))
                error 'new header has different sublabel(s) from current header'
            end
            
            % that's it for 'all' flag
            if strcmp(flag, 'all'), return, end
            
            % check number of points
            switch flag
                case 'new'
                    n_check = H.n + length(ind);
                case 'remove'
                    n_check = H.n - length(ind);
                case 'chg&new'
                    n_check = H.n + length(ind{2});
                case 'chg&rm'
                    n_check = H.n - length(ind{2});
                otherwise
                    n_check = H.n;
            end
            if new_head.n ~= n_check
                error 'wrong number of elements in new header'
            end
            
            % check values
            if H.is_measure
                switch flag
                    case 'all'
                        % nothing more needs to be checked
                    case 'chg'
                        % change is possible only if all values were
                        % changed
                        ok = (new_head.start == H.start && new_head.scale == H.scale) || (length(ind) == H.n);
                    case {'chg&new', 'new', 'chg_data'}
                        ok = (new_head.start == H.start) && (new_head.scale == H.scale);
                    case {'remove', 'chg&rm'}
                        if strcmp(flag, 'remove'), idx_rm = ind; else idx_rm = ind{2}; end
                        ok = (new_head.scale == H.scale);
                        if ok
                            if idx_rm(1) == 1 && idx_rm(end) == H.n
                                nrm_start = find(diff(idx_rm) ~= 1);
                                if isempty(nrm_start)
                                    % removed all!
                                    ok = (new_head.start == H.start);
                                elseif isscalar(nrm_start)
                                    ok = abs(new_head.start - (H.start + nrm_start*H.scale)) < abs(H.scale)/1000;
                                else
                                    ok = false;
                                end
                            elseif idx_rm(1) == 1
                                ok = abs(new_head.start - (H.start + length(idx_rm)*H.scale)) < abs(H.scale)/1000 && all(diff(idx_rm) == 1);
                            elseif idx_rm(end) == H.n
                                ok = (new_head.start == H.start) && all(diff(idx_rm) == 1);
                            else
                                ok = false;
                            end
                        end
                    case 'perm'
                        ok = false;
                    otherwise
                        error('invalid flag ''%s'' for header update', flag)
                end
            else
                switch flag
                    case 'all'
                        % nothing more needs to be checked
                        [idx_h, idx_n] = deal([]);
                    case 'chg'
                        unchg = true(1, H.n);
                        unchg(ind) = false;
                        [idx_h, idx_n] = deal(find(unchg)); % (much) faster than setdiff
                    case {'new', 'chg_data', 'sub_data'}
                        [idx_h, idx_n] = deal(1:H.n);
                    case 'chg&new'
                        unchg = true(1, H.n);
                        unchg(ind{1}) = false;
                        [idx_h, idx_n] = deal(find(unchg));
                    case 'chg&rm'
                        unchg = true(1, H.n);
                        unchg([ind{:}]) = false;
                        idx_h = find(unchg);
                        unchg = true(1, new_head.n);
                        unchg([ind{1}]) = false;
                        idx_n = find(unchg);
                    case 'remove'
                        unchg = true(1, H.n);
                        unchg(ind) = false;
                        idx_h = find(unchg);
                        idx_n = 1:new_head.n;
                    case 'perm'
                        idx_h = ind;
                        idx_n = 1:H.n;
                    otherwise
                        error('invalid flag ''%s'' for header update', flag)
                end
                ok = isequal(new_head.values(idx_n, 1:H.n_column), H.values(idx_h, :));
            end
            if ~ok, error 'New header values are not consistent with the operation', end
        end
        function new_head = add_label(H, label, value)
            % xplr.DimensionLabel object
            if isa(label, 'xplr.DimensionLabel')
                name = label.label;
                if nargin < 3, value = label.default_val; end
            else
                name = label;
                if iscell(value), x = value{1}; else x = value; end % there must be a value to infer type correctly
                type = xplr.DimensionLabel.infer_type(x);
                label = xplr.DimensionLabel(name, type);
            end
            if any(strcmp(name, {H.sub_labels.label}))
                error('a sub-label ''%s'' already exists', name)
            end
            % add it to the list of labels
            idx = length(H.sub_labels) + 1;
            new_head = copy(H);
            new_head.sub_labels(idx) = label;
            if iscell(value)
                new_head.values(:, idx) = value;
            else
                if size(new_head.values, 1) == 0
                    new_head.values = cell(0, idx);
                else
                    [new_head.values{:,idx}] = deal(value);
                end
            end
        end
        function new_values = track_values(H, indices)
            n_ind = length(indices);
            new_values = cell(n_ind, H.n_column);
            for i=1:n_ind
                idx = indices{i};
                v = H.values(idx, :);
                if isscalar(idx)
                    % just copy values
                    new_values(i, 1:H.n_column) = v;
                else
                    % several indices together
                    for j=1:H.n_column
                        switch H.sub_labels(j).label
                            case 'ViewColor'
                                % average colors
                                new_values{i, j} = mean(cat(1, v{:,j}), 1);
                            case 'Name'
                            otherwise
                                % default behavior: compute union of values
                                switch H.sub_labels(j).type
                                    case {'numeric', 'logical'}
                                        new_values{i, j} = unique([v{:,j}], 'stable');
                                    case 'char'
                                        new_values{i, j} = unique(v(:,j)', 'stable');
                                end
                        end
                    end
                end
            end
        end
    end
    
end
