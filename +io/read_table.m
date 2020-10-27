function data = read_table(file)
% function data = read_table(file)

if nargin<1
    file = brick.getfile('*.csv', 'Select CSV file');
end

file = cellstr(file);
n_file = length(file);

% Data name
if n_file > 1
    file_base = brick.map(@(f)brick.fileparts(f, 'base'), file);
    name = char(file_base);
    all_same = ~any(diff(name));
    name = name(1,:); name(~all_same) = '*';
else
    name = brick.fileparts(file{1}, 'base');
end

% Fix import options:
% if numeric values are surrounded by quotes, read them as string!
f = file{1};
opts = detectImportOptions(f);
fid = fopen(f, 'r');
for k = 1:opts.DataLines(1)
    line = fgetl(fid);
end
fclose(fid);
line = strsplit(line, opts.Delimiter);
for k = find(strcmp(opts.VariableTypes, 'double'))
    c = line{k}(1);
    if c == '''' || c == '"'
        opts.VariableTypes{k} = 'string';
    end
end

% Read table for first table
a = readtable(f, opts);

p = a.Properties;
column_names = p.VariableNames;
[n_row, n_col] = size(a);

% Detect dates and times where Matlab maybe didn't
v_types = opts.VariableTypes;
column_convert_datetime = false(1, n_col);
for i = find(strcmp(v_types, 'char'))
    x = a{1, i};
    % remove everything after sign '+' and try convert to datetime
    x = regexprep(x, '+.*', '');
    try
        datetime(x);
        % succeeded in converting to datetime! -> convert the full column
        column_convert_datetime(i) = true;
        a.(i) = datetime(regexprep(a{:, i}, '+.*', ''));
    catch
    end
end

% So far we have a 2D table, dimensions are:
% - columns = "variables"
% - rows = "samples"
% Detect which column bears the samples header info (strings, dates, ids,
% repetition number etc, as opposed to the samples numerical data), and
% whether data is sampled in a regular manner (i.e. 1 data sample per
% header value) so that it can be reshaped as a multi-dimensional array.

% Scan columns to detect additional dimensions that can be unfolded.
% In particular, detect "inside repetition" and "outside repetition"; for
% example if one column as values 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, etc.,
% then "internal_repeat" is 2 and "n_values" is 3.
[values, header_spec] = deal(cell(1, n_col));
[internal_repeat, n_values] = deal(zeros(1, n_col)); 
[is_header, is_data] = deal(false(1, n_col));
for k = 1:n_col
    [values_k, ~, index] = unique(a(:, k), 'stable');
    values_k = values_k.Variables; % values of column i
    nk = length(values_k);
    if nk == 1
        % same value everywhere! we are not interested in this column
        internal_repeat(k) = n_row;
        is_header(k) = false;
        is_data(k) = isnumeric(values_k);
        continue
    elseif nk == n_row
        % no repetition
    else
        % some repetitions: detect systematic internal/external repetition
        % if any
        diff_with_next = [logical(diff(index)); true];
        internal_rep = diff(find(diff_with_next));
        if any(diff(internal_rep))
            internal_rep = 1;
        else
            internal_rep = internal_rep(1);
        end
        if mod(n_row, nk * internal_rep)
            external_rep = 1;
        else
            external_rep = n_row / (nk * internal_rep);
        end
        if (internal_rep ~= 1 || external_rep ~= 1) ...
                && isequal(index, ...
                repmat(kron((1:nk)', ones(internal_rep, 1)), [external_rep 1]))
            internal_repeat(k) = internal_rep;
            n_values(k) = nk;
        end
    end
    
    % detect measure
    isnum = isnumeric(values_k) || isdatetime(values_k) || isduration(values_k);
    if isnum
        d = diff(values_k);
        equally_spaced = all(abs(diff(d)) < d(1)*1e-6);
    else
        equally_spaced = false;
    end
    if equally_spaced
        % header spec = unit, n, scale, start
        header_spec{k} = {'', nk, d(1), values_k(1)};
    else
        header_spec{k} = {values_k};
    end

    % detect whether values are likely to be header info
    is_header(k) = internal_repeat(k) || ~isnum || equally_spaced;
    is_data(k) = ~is_header(k);
end

% Header for rows/samples
% -> Attempt to exploit repetition to unfold multidimensional array
sub_dims = find(is_header & internal_repeat);
if ~isempty(sub_dims)
    % Detected sub-dimensions
    cycle_length = internal_repeat(sub_dims) .* n_values(sub_dims);
    [cycle_length, dim_order] = sort(cycle_length);
    sub_dims = sub_dims(dim_order);
    % the internal repetitions should be predicted by the cycles of other
    % dimensions
    ok = (internal_repeat(sub_dims(1)) == 1) ...
        && all(internal_repeat(sub_dims(2:end)) == cycle_length(1:end-1)) ...
        && cycle_length(end) == n_row;
    if ~ok
        error 'Failure in detecting sub-dimensions'
    end
    % there should be no other header column
    if any(is_header & ~internal_repeat)
        error 'Detected sub-dimensions but other column(s) seem to bear header information'
    end
    rows_header = cell(1, length(sub_dims));
    for i = 1:length(sub_dims)
        k = sub_dims(i); % column number
        rows_header{i} = xplr.Header(column_names{k}, header_spec{k}{:});
    end
    rows_header = [rows_header{:}];
else
    % Did not detect sub-dimensions
    header_dim = find(is_header);
    if isempty(header_dim)
        % did not detect headers, use a mere enumeration
        rows_header = xplr.Header(name, n_rows);
    elseif isscalar(header_dim)
        % use specification (can be e.g. for a measure header)
        rows_header = xplr.Header(column_names{header_dim}, header_spec{header_dim}{:});
    else
        % table header
        rows_header = xplr.Header(column_names(header_dim), table2cell(a(:, header_dim)));
    end
end

% Header for columns and data
if ~any(is_data)
    error 'no data'
elseif sum(is_data) == 1
    % one-dimension data, no header for columns
    columns_header = xplr.Header.empty(1, 0);
else
    columns_header = xplr.Header('Variables', column_names(is_data));
end

% Header for files
if n_file == 1
    files_header = xplr.Header.empty(1, 0);
else
    files_header = xplr.Header('File', file_base);
end

% Data -> read from all files only now!
dat = zeros([n_row, sum(is_data), n_file]);
for k = 1:n_file
    if k == 1
        ak = a;
    else
        % Read table from next files
        ak = readtable(file{k}, opts);

        % Checkups
        pk = ak.Properties;
        if size(ak, 1) ~= n_row
            error(['File ''' file_base{1} ''' has ' n_row ' samples,' ...
                ' but file ''' file_base{k} ''' has ' num2str(size(ak, 1)) ' samples.'])
        elseif size(ak, 2) ~= n_col
            error(['File ''' file_base{1} ''' has ' n_col ' columns,' ...
                ' but file ''' file_base{k} ''' has ' num2str(length(pk.VariableNames)) ' columns.'])
        elseif ~isequal(pk.VariableNames, column_names)
            problem = false(1, n_col);
            for i = 1:n_col
                problem(i) = ~strcmp(pk.VariableNames{i}, column_names{i});
            end
            error(sprintf( ...
                ['File ''' file_base{1} ''' and ''' file_base{k} ''' do not have the ' ...
                'same column names for column(s) ' fn_strcat(find(problem), ', ') '.\n' ...
                file_base{1} ': ' fn_strcat(column_names(problem), ', ') '.\n' ...
                file_base{k} ': ' fn_strcat(pk.VariableNames(problem), ', ') '.']))
        end

        % Detect dates and times where Matlab maybe didn't
        for i = find(column_convert_datetime)
            ak.(i) = datetime(regexprep(ak{:, i}, '+.*', ''));
        end

        % Header identity
        assert(isequal(ak(:, is_header), a(:, is_header)))
    end
    dat(:, :, k) = table2array(ak(:, is_data));
end

% Replace 0 values by NaN
z = (dat == 0);
if any(z(:))
    answer = questdlg('Replace zero values by NaNs?', 'XPLOR', ...
        'Yes', 'No', 'Yes');
    if strcmp(answer, 'Yes')
        dat(z) = NaN;
    end
end

% Unfold data (i.e. reshape if there are sub-dimensions)
if any(sub_dims)
    dat = reshape(dat, [n_values(sub_dims), sum(is_data), n_file]);
end


% Separate date and time !!
header = [rows_header columns_header files_header];
nd = length(header);
for k = 1:nd
    head_k = header(k);
    if head_k.is_datetime
        day_scale = days(head_k.scale); % double
        day_span = head_k.n * day_scale;
        if day_span > 4 && day_scale < .2 && (mod(1, day_scale) == 0)
            % separate data between days and times!!
            day_start = dateshift(head_k.start, 'start', 'day');
            day_header = xplr.Header('Day', '', ceil(day_span), days(1), day_start);
            time_header = xplr.Header('Time', '', 1 / day_scale, head_k.scale, day_start);
            % new header
            header = [header(1:k-1) time_header day_header header(k+1:end)];
            % new data
            sz = size(dat);
            sz2 = [sz(1:k-1), time_header.n day_header.n sz(k+1:end)];
            sz_ = [prod(sz(1:k-1)), head_k.n prod(sz(k+1:end))];
            sz2_ = [prod(sz(1:k-1)), time_header.n*day_header.n prod(sz(k+1:end))];
            dat2_ = zeros(sz2_, 'like', dat);
            offset = days(head_k.start - day_start) / day_scale;
            dat_ = reshape(dat, sz_);
            dat2_(:, offset + (1:head_k.n), :) = dat_(:, :, :);
            dat = reshape(dat2_, sz2);
            break
        end
    end
end

% Structured data
data = xplr.XData(dat, header, name);

