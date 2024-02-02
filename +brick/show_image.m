function show_image(img, background_color)
% Display an image with alpha channel and with first 2 dimensions being x-y
%---
% function show_image(img [,color|'checker'])
%---
% Input:
% - img     nx*ny*nc array, with nc being 1, 3 or 4
% - color   a color flag ('k', 'r', etc.), a 1x3 color value, or 'checker':
%           the background color to use (default: 'checker')
%
% See also: brick.white2alpha

% Convert uint to double
switch class(img)
    case 'uint8'
        img = single(img) / 255;
    case 'uint16'
        img = single(img) / 65535;
    case {'single' 'double' 'logical'}
        % ok
    otherwise
        error 'type not handled'
end

% Add background for image with alpha channel
[nx, ny, nc] = size(img);
if nc == 4
    alpha = img(:,:,4);
    img = img(:,:,1:3);
    if nargin < 2, background_color = 'checker'; end
    if strcmp(background_color, 'checker')
        step = round(mean([nx ny])/min([nx ny 15]));
        xstripes = mod(floor((0:nx-1)'/step),2);
        ystripes = mod(floor((0:ny-1)/step),2);
        background = bsxfun(@xor,xstripes,ystripes);
        background = .45 + .1*background;
    else
        background_color = brick.colorbyname(background_color);
        background = brick.third(background_color);
    end
    img = brick.add(brick.mult(background, 1-alpha), brick.mult(img, alpha));
end

% Display
imagesc(permute(img, [2 1 3]))
axis image
colormap gray