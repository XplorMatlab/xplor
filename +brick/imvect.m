function x = imvect(x,mask,varargin)
%IMVECT Convert an image to a vector of pixels inside a mask, and vice-versa
%---
% function x = imvect(x,mask[,outputtype][,outsidevalue])
% function x = imvect(x[,'vector'])
% function k = imvect(ij,mask[,outsidevalue])
% function ij = imvect(k,mask[,outsidevalue])
%---
% switch between "image" and "vector" representation of the pixels in an
% image
%
% Input:
% - x       array of size (nx,ny,nt,...) or (np,nt,...)
% - mask    logical array of size (nx,ny), such that sum(mask(:))==np,
%           or dimensions [nx ny] to create the mask true(nx,ny) 
%           mask can also be a cell array of logical arrays: then operation
%           applies to succesive dimensions
% - outputtype      'vector' or 'image': default behavior toggles
%                   representation; by setting outputtype, x is unchanged if
%                   is already has the desired representation
%                   'maskimage': if input is the image, output will be the
%                   masked image
% - outsidevalue    value to set outside the mask in the image [default=0]
% 
% Output:
% - x       array size became (np,nt,...) or (nx,ny,nt,...), respectively
%
% See also brick.indices, brick.maskavg, brick.maskselect

% Thomas Deneux
% Copyright 2011-2017

if nargin==0, help brick.imvect, return, end

% input
if isvector(x), x = x(:); end
if nargin<2 || strcmp(mask,'vector'), mask = [size(x,1),size(x,2)]; end
if ~iscell(mask), mask = {mask}; end
[nx ny np] = deal(zeros(1,length(mask)));
for i=1:length(mask)
    maski = mask{i};
    if ~islogical(maski) && isvector(maski) && length(maski)==2
        nx(i) = maski(1); ny(i) = maski(2);
        np(i) = nx(i)*ny(i);
        mask{i} = [];
    else
        if ~islogical(maski)
            error('mask(s) must be a logical array or a lenth 2 vector')
        end
        [nx(i) ny(i)] = size(maski);
        np(i) = sum(maski(:));
    end
end
outsidevalue = 0;
outputtype = 'toggle';
for k=1:nargin-2
    a = varargin{k};
    if isnumeric(a)
        outsidevalue = a;
    else
        if ~brick.ismemberstr(a,{'image' 'maskimage' 'vector' 'toggle'}), error('invalid flag ''%s''',a), end
        outputtype = a;
    end
end

% special: indices
if isscalar(x) || (isvector(x) && length(x)==2)
    if ~isscalar(mask), error 'argument', end % note that nx, ny and np are scalars as well
    mask = mask{1};
    if isscalar(x)                              % vector index to image index
        k = x;
        test = zeros(np,1); test(k) = 1;
        test = brick.imvect(test,mask);
        [i j] = find(test);
        x = [i j];
    elseif isvector(x) && length(x)==2          % image index to vector index
        ij = x;
        i = ij(1); j = ij(2);
        test = false(nx,ny); test(i,j) = true;
        test = brick.imvect(test,mask);
        k = find(test);
        x = k;
    end
    return
end

% output type
s = size(x);
if all(s(1:2)==[nx(1) ny(1)])                   % image to vector
    if strcmp(outputtype,'image')
        return % (nothing to do)
    elseif strcmp(outputtype,'maskimage')
        % special
        x = brick.imvect(brick.imvect(x,mask,'vector'),mask,'image',outsidevalue);
        return
    end
    outputtype = 'vector';
elseif s(1)==np(1)                              % vector to image
    if strcmp(outputtype,'vector'), return, end % (nothing to do)
    outputtype = 'image';
else
    error('dimensions of x and mask do not fit')
end

% loop on masks
for i=1:length(mask)
    s = size(x);
    maski = mask{i};
    switch outputtype
        case 'vector'
            s1 = s(1:(i-1)); s2 = s(i+2:end);
            x = reshape(x,[prod(s1) nx(i)*ny(i) prod(s2)]);
            if ~isempty(maski), x = x(:,maski,:); end
            x = reshape(x,[s1 np s2 1]);
        case 'image'
            s1 = s(1:2*(i-1)); s2 = s(2*i:end);
            x = reshape(x,[prod(s1) np prod(s2)]);
            if ~isempty(maski)
                xmem = x;
                if outsidevalue==0
                    x = zeros([prod(s1) nx(i)*ny(i) prod(s2)],class(x)); %#ok<*ZEROLIKE>
                elseif outsidevalue==1
                    x = ones([prod(s1) nx(i)*ny(i) prod(s2)],class(x));
                elseif isnan(outsidevalue)
                    x = nan([prod(s1) nx(i)*ny(i) prod(s2)],class(x));
                else
                    x = outsidevalue*ones([prod(s1) nx(i)*ny(i) prod(s2)],class(x));
                end
                x(:,maski,:) = xmem;
            end
            x = reshape(x,[s1 nx(i) ny(i) s2]);
    end
end
