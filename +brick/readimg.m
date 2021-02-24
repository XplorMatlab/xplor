function [a, alpha] = readimg(fname,varargin)
% function [a [,alpha]] = readimg(fname[,'nopermute'][,'double'])
%---
% read image using imread, and handles additional features:
% - converts to double if 'double' option is set
% - detects if color or gray-scale images (in the last case, use a 2D array per image)
% - can read a stack of images (returns 3D array)
%
% images are read according to x-y convention, use 'nopermute' flag to use
% Matlab y-x convention

% Thomas Deneux
% Copyright 2004-2017

% Input
if nargin<1 || isempty(fname)
    fname = brick.getfile;
end
[nopermute, dofloat] = brick.flags('nopermute', 'double', varargin);
dopermute = ~nopermute;
fname = cellstr(fname);
nimages = length(fname);

% first image
[a, alpha] = readoneframe(fname{1}, dopermute, true, false);
firstchannelonly = (size(a,3) == 1);

% multi-gif?
if size(a,4) > 1 && nimages > 1
    error('cannot read multiple GIF files')
end

% multi-tiff?
if strfind(lower(brick.fileparts(fname{1},'ext')),'.tif')
    nframes = length(imfinfo(fname{1}));
    if nimages>1 && nframes>1, error 'cannot handle multiple tif that themselves have multiple frames', end
else
    nframes = 1;
end
        
% stack
if nimages*nframes>1
	brick.progress('reading frame',nimages*nframes)
    if firstchannelonly
        a(1,1,nimages*nframes) = 0;
    else
        a(1,1,1,nimages*nframes) = 0;
    end
    if ~isempty(alpha)
        alpha(1,1,1,nimages*nframes) = 0;
    end
    for i=2:nimages*nframes
        brick.progress(i)
        if nframes>1
            [b beta] = readoneframe(fname{1},dopermute,false,firstchannelonly,i);
        else
            [b beta] = readoneframe(fname{i},dopermute,false,firstchannelonly);
        end        
        if firstchannelonly
            a(:,:,i) = b(:,:,1);
        else
            a(:,:,:,i) = b;
        end
        if ~isempty(alpha)
            alpha(:,:,i) = beta;
        end
    end
    brick.progress('end')
end

% make float-encoded color image btw 0 and 1
if dofloat || ismember(class(a), {'single', 'double'})
    switch class(a)
        case {'single' 'double'}
            nbyte = ceil(log2(max(a(:)))/8);
        case 'uint8'
            nbyte = 1;
        case 'uint16'
            nbyte = 2;
        otherwise
            if brick.dodebug, disp 'please help me', keyboard, end
    end
    a = brick.float(a);
    alpha = brick.float(alpha);
    switch nbyte
        case 0
            % max(a(:)) is 1, this is fine
        case 1
            a = a/255;
            alpha = alpha/255;
        case 2
            a = a/65535;
            alpha = alpha/65535;
        otherwise
            if brick.dodebug, disp 'please help me', keyboard, end
    end
end

% combine alpha with image?
if nargout < 2 && ~isempty(alpha)
    if firstchannelonly
        a = repmat(a, [1 1 3]);
    end
    a = cat(3, a, alpha);
end
    

%---
function [a alpha] = readoneframe(f, dopermute, docheckgrayscale, firstchannelonly, i)

if nargin<5
    [a cmap alpha] = imread(f); 
else
    [a cmap alpha] = imread(f, i);
end
if dopermute
    a = permute(a,[2 1 3 4]); % Matlab (y,x) convention -> convention (x,y)
    if ~isempty(alpha)
        alpha = permute(alpha, [2 1 3]);
    end
end

% we can have multiple images for gif files
multiimage = (size(a,4) > 1);

if size(a,3)==3
    firstchannelonly = firstchannelonly || (docheckgrayscale &&  ~any(any(any(diff(a,1,3)))));
    if firstchannelonly
        % detected grayscale image saved as a color image: keep only one
        % channel since the three of them are identical
        a = a(:,:,1,:);
    end
elseif ~isempty(cmap)
    % apply colormap
    firstchannelonly = firstchannelonly || (docheckgrayscale &&  ~any(any(diff(cmap,1,2))));
    if firstchannelonly
        % grayscale colormap!
        a = reshape(cmap(a(:)+1,1), size(a));
    else
        a = reshape(cmap(a(:)+1,:),[size(a) 3]);
        if multiimage
            % a is nx*ny*1*nfr*3 -> make it nx*ny*3*nfr
            a = permute(a,[1 2 5 4 3]);
        end
    end
end

