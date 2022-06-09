function indices = color2index(a, cmap, hsv_weight)
% Convert color image to indices in the colormap
% function indices = color2index(a, cmap[, hsv_weight])
%---
% Input:
% - a       color image
% - cmap    color map: Nx3 array or string
% - hsv_weight  value between 0 and 1; work in RGB space if 0, in HSV space
%           if 1, inbetween for intermediary values
%
% Output:
% - indices image of indices in color map pixel colors in a are the closest
%           to
%
% If color map is specified as values, the difference between color values
% is computed. If color map is specified as a string, very particular
% calculation related to known color maps are undertaken.
%
% See also brick.color2bw, brick.index2color

% Input -> make sure color are values between 0 and 1
switch class(a)
    case {'float' 'double'}
        assert(all(a(:)>=0 & a(:)<=1), 'Color image values must be between 0 and 1.')
    case 'uint8'
        a = double(a) / 255;
    otherwise
        error('input image of class %s not handled yet', class(a))
end
switch class(cmap)
    case {'float' 'double'}
        assert(all(cmap(:)>=0 & cmap(:)<=1), 'Color values in color map must be between 0 and 1.')
    case 'uint8'
        cmap = double(cmap) / 255;
    case {'char' 'string'}
        cmap = char(cmap);
        if ~ismember(cmap, {'mapgeog'})
            error('colormap "%s" not handled yet', cmap)
        end
    otherwise
        error('input image of class %s not handled yet', class(a))
end
if nargin<3
    hsv_weight = 0;
end


% Very different processing depending on whether color map is numerical or
% a string
if isnumeric(cmap)
    % Add HSV values in addition to RGB values
    if hsv_weight > 0
        a2 = rgb2hsv(a);
        cmap2 = rgb2hsv(cmap);
        if hsv_weight == 1
            [a, cmap] = deal(a2, cmap2);
        else
            a = (1-hsv_weight) * a + hsv_weight * a2;
            cmap = (1-hsv_weight) * cmap + hsv_weight * cmap2;
        end
    end

    % Distances to colormap
    diff = a - permute(cmap,[3 4 2 1]);
    dist = sum(diff.^2, 3);
    [~, indices] = min(dist, [], 4);

else
    % 
    
end

