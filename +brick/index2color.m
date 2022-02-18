function a = index2color(indices, cmap)
% Convert indices to color image (use brick.clip if any scaling is needed)
% function a = index2color(indices, cmap)
%---
% Input:
% - indices image of indices in the color map
% - cmap    color map: Nx3 array or string
%
% Output:
% - a       color image
%
%
% See also brick.clip, brick.color2index

% Input -> make sure color are values between 0 and 1
switch class(indices)
    case {'float' 'double'}
    case 'uint8'
        indices = double(indices) + 1;
    otherwise
        error('input image of class %s not handled yet', class(a))
end
switch class(cmap)
    case {'float' 'double'}
        assert(all(cmap(:)>=0 & cmap(:)<=1), 'Color values in color map must be between 0 and 1.')
    case 'uint8'
        cmap = double(cmap) / 255;
    case {'char' 'string'}
        try
            cmap = eval(cmap);
        catch
            cmap = eval(strcat("colormaps.", cmap));
        end
    otherwise
        error('input image of class %s not handled yet', class(a))
end

% Color image
a = reshape(cmap(indices(:), :), [size(indices) 3]);