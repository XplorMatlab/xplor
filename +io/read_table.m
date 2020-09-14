function data = read_table(file)
% function data = read_table(file)

file = cellstr(file);
n_file = length(file);

% Detect type
data = cell(1, n_file);
for i = 1:n_file
    data{i} = read_one_file(file{i});
end
data = [data{:}];

%---
function data = read_one_file(f)
%%

% Data name
name = fn_fileparts(f, 'base');

% Fix import options:
% if numeric values are surrounded by quotes, read them as string!
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

% Read table
a = readtable(f, opts);

p = a.Properties;
column_names = p.VariableNames;
[n_row, n_col] = size(a);

% Detect dates and times where Matlab maybe didn't
v_types = opts.VariableTypes;
for i = find(strcmp(v_types, 'char'))
    x = a{1, i};
    % remove everything after sign '+' and try convert to datetime
    x = regexprep(x, '+.*', '');
    try
        datetime(x);
        % succeeded in converting to datetime! -> convert the full column
        opts.VariableTypes{i} = 'datetime';
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
% -> Attempt to exploit repetition to build-up multidimensional array
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
        rows_header = xplr.Header(column_names(header_dim), header_spec{header_dim}{:});
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
    columns_header = xplr.Header('data', column_names(is_data));
end

% Data
dat = table2array(a(:, is_data));
if any(sub_dims)
    dat = reshape(dat, [n_values(sub_dims), sum(is_data)]);
end

% Structured data
data = xplr.XData(dat, [rows_header columns_header], name);

