% flow dataset
[~, ~, ~, dat] = flow;
dat = permute(dat, [1, 3, 2]);
s = size(dat);

% display
xplor(dat, ...
    'header', {{'x', 'px'}, {'y', 'px'}, {'axis', num2cell(['a':'y' 'A':'Y'])}}, ...
    'name', 'Flow Data');

