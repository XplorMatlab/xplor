function [b2, beta2] = combine_images(a, b, out_type)
% Superimpose image a on top of b (both with transparency)
%---
% function b = combine_images(a, b)
% function [b, alpha] = combine_images(a, b [,'uint8|double'])
%---
% Input:
% - a, b    Images with transparency: a 2D array (1-channel, grayscale
%           image, w), 2-channel (wa) or 3-channel (rgb) or 4-channel
%           (rgba)3D array, or a cell array {img, alpha}.
% - 'uint8|double'  Whether output should be uint8 (between 0 and 255) or
%           double (between 0 and 1). If not specified, using the same
%           format as for input b.
%
% Output:
% - b2      New image with a superimposed onto the input b. Normally a
%           4-channel rgba image, except if the alpha channel is requested
%           to be output to a second output 'beta', in which case b will be
%           a 2D array if both a and b are grayscale, or a 3-channel rgb
%           image otherwise.
% - beta2   Alpha channel of the output.

% Input
% (decompose alpha channel)
[a, alpha] = decompose_alpha(a);
[b, beta] = decompose_alpha(b);
% (output type)
if nargin < 3
    out_type = class(b);
end
% (convert to double)
a = image2float(a);
alpha = image2float(alpha);
b = image2float(b);
beta = image2float(beta);

% Superimpose a onto b

%   superimposing b on x = beta * b + (1-beta) * x
%   superimposing a on the result = alpha * a + (1-alpha) * (b * beta + x * (1-beta))
%                                 = (alpha*a + beta*(1-alpha)*b) + (1-alpha)*(1-beta) * x
%                                 = b2 * beta2 + (1-beta2) * x

% beta2 = 1 - (1-alpha).*(1-beta);
beta2 = alpha + beta .* (1-alpha);
b2 = alpha .* a + beta .* (1-alpha) .* b;
b2 = b2 ./ beta2;
b2(isnan(b2)) = 0;

% Output
[nx, ny, nc] = size(b2);
if isscalar(beta2)
    % happens when both alpha and beta were scalars!
    beta2 = beta2 * ones(nx, ny);
end
switch out_type
    case 'uint8'
        b2 = uint8(b2 * 255);
        beta2 = uint8(beta2 * 255);
    case 'uint16'
        b2 = uint8(b2 * 65535);
        beta2 = uint8(beta2 * 65535);
end
if nargout < 2
    if nc == 1
        b2 = repmat(b2, [1 1 3]);
    end
    b2 = cat(3, b2, beta2);
end

%---
function [a, alpha] = decompose_alpha(a)

if iscell(a)
    [a, alpha] = deal(a{:});
elseif ismember(size(a,3), [2 4])
    [a, alpha] = deal(a(:,:,1:end-1), a(:,:,end));
elseif ismember(size(a,3), [1 3])
    % no transparency!
    alpha = 1;
end

%--- 
function a = image2float(a)

switch class(a)
    case {'single', 'double'}
    case 'uint8'
        a = single(a) / 255;
    case 'uint16'
        a = single(a) / 65535;
    otherwise
        error('wrong input class')
end

