function colors = colors(name, k)
% function colors = colors(colorname)
% function colors = colors(set_name, k)
%---
% Returns specified color, or at indices k from specific map. 
% Available maps are 'link_key'. Some maps are infinite: they start with a
% set of specific pre-determined values, and continue with random values.

% Thomas Deneux
% Copyright 2008-2012

% named color or index in set?
if nargin == 1
    tokens = fn_strcut(name,'.');
    colors = [];
    switch tokens{1}
        case 'gui'
            switch tokens{2}
                case 'controls'
                    switch tokens{3}
                        case 'dataname'
                            colors = .6;
                        case 'filter'
                            colors = .85;
                        case 'item'
                            colors = 1;
                    end
            end
    end
    if isempty(colors), error('unknown color name ''%s''', name), end
    if isscalar(colors), colors = colors([1, 1, 1]); end
    return
end

% color map
% tokens = regexp(set_name,'^([^\d]*)(\d*)$','tokens');
% set_name = tokens{1}{1};
% ncol = str2double(tokens{1}{2});
infinite_map = false;
switch name
    case 'link_key'
        sat_range = .3;
        lum_range = 1;
        c_map = squeeze(hsv2rgb([4, 6, 2, 1, 3]/6, ones(1, 5)*sat_range,ones(1, 5)*lum_range));
        infinite_map = true;
        % index might be zero! add one to make indexing start from 1
        c_map = [[1, 1, 1]*.94; c_map];
        k = k + 1;
    otherwise
        error('unknown color set name ''%s''', name)
end

if infinite_map
    n_map = size(c_map, 1);
    
    ok_map = (k <= n_map);
    if all(ok_map)
        colors = c_map(k, :);
    else
        colors = zeros(length(k), 3);
        colors(ok_map, :) = c_map(k(ok_map), :);
        k_random = k(~ok_map) - n_map;
        c_map_random = randomcolors(max(k_random), sat_range, lum_range);
        colors(~ok_map, :) = c_map_random(k_random, :);
    end
end
