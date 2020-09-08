function data = read_table(file)
% function data = read_table(file)

file = cellstr(file);
n_file = length(file);

data = struct('header', cell(1, n_file), ...
    'data', cell(1, n_file));
for i = 1:n_file
    [data(i).header, data(i).data] = read_one_file(file{i});
end

%---
function [header, data] = read_one_file(f)

a = read_table(f);
%%
p = a.Properties;
column_names = p.VariableNames;
[n_row, n_col] = size(a);

%%
% Scan columns to detect header dimensions

% detect "inside repetition", "outside repetition"; for example if one
% column as values 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, etc., then "inside
% repetition" is 2 and "outside repetition" is 3
internal_repeat = zeros(1, n_col); 
for i = 1:n_col
    [values, i_first, index] = unique(a(:, i));
    if length(values) == 1
        % same value everywhere!
        internal_repeat(i) = n_row;
    elseif length(values) == n_row
        % no repetition
        internal_repeat(i) = 0;
    else
        % some repetitions: detect systematic internal repetition if any
        diff_with_next = [any(diff(index)); True];
        dist = diff(find(diff_with_next)); % distance with next diff
        rep = min(dist);
        if rep > 1 && ~any(mod(dist, rep))
            internal_repeat(i) = rep;
        end
    end
end
