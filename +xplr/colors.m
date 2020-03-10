function colors = colors(name,k)
% function colors = colors(colorname)
% function colors = colors(setname,k)
%---
% Returns specified color, or at indices k from specific map. 
% Available maps are 'linkkey'. Some maps are infinite: they start with a
% set of specific pre-determined values, and continue with random values.

% Thomas Deneux
% Copyright 2008-2012

% named color or index in set?
if nargin==1
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
    if isempty(colors), error('unknown color name ''%s''',name), end
    if isscalar(colors), colors = colors([1 1 1]); end
    return
end

% color map
% tokens = regexp(setname,'^([^\d]*)(\d*)$','tokens');
% setname = tokens{1}{1};
% ncol = str2double(tokens{1}{2});
infinitemap = false;
switch name
    case 'linkkey'
        satrange = .3;
        lumrange = 1;
        cmap = squeeze(hsv2rgb([4 6 2 1 3]/6,ones(1,5)*satrange,ones(1,5)*lumrange));
        infinitemap = true;
        % index might be zero! add one to make indexing start from 1
        cmap = [[1 1 1]*.94; cmap];
        k = k+1;
    otherwise
        error('unknown color set name ''%s''',name)
end

if infinitemap
    nmap = size(cmap,1);
    
    okmap = (k<=nmap);
    if all(okmap)
        colors = cmap(k,:);
    else
        colors = zeros(length(k),3);
        colors(okmap,:) = cmap(k(okmap),:);
        krandom = k(~okmap)-nmap;
        cmaprandom = randomcolors(max(krandom),satrange,lumrange);
        colors(~okmap,:) = cmaprandom(krandom,:);
    end
end
