function a = readmovie(filename,varargin)
% function a = readmovie(filename,frames[,'nodisplay'][,'nopermute'],['bin',[xbin tbin]])
%---
% read an avi file and stores it into a 2D+time array 
% (2D+time+channel if color movie)
%
% See also brick.savemovie

% Thomas Deneux
% Copyright 2004-2017

if nargin<1 || isempty(filename)
    filename = brick.getfile('*.avi');
end
if ~exist(filename,'file')
    error('file ''%s'' does not exist',filename)
end
frames = {};
dodisplay = true; dopermute = true; bin = [1 1];
k = 1;
while k <= length(varargin)
    a = varargin{k};
    k = k+1;
    if isnumeric(a)
        frames = a;
    elseif ischar(a)
        switch a
            case 'nodisplay'
                dodisplay = false;
            case 'nopermute'
                dopermute = false;
            case 'bin'
                bin = varargin{k};
                k = k+1;
            otherwise
                error argument
        end
    else
        error argument
    end
end

if dodisplay, disp 'reading', end
% try
    % recent Matlab version
    if ~isempty(frames), frames = {frames([1 end])}; end
    v = VideoReader(filename);
    if all(bin == 1)
        a = read(v,frames{:},'native');
    else
        xbin = bin(1);
        tbin = bin(2);
        w = floor(v.Width / xbin);
        h = floor(v.Height / xbin);
        nframe = floor(v.NumFrames / tbin);
        n_per_block = round(100 / tbin);
        nblock = floor(nframe / n_per_block);
        nframe = nblock * n_per_block;
        img = v.read(1,'native');
        nc = size(img,3);
        a = zeros(w, h, nc, nframe, class(img));
        brick.progress('reading block', nblock)
        for kblock = 1:nblock
            brick.progress(kblock)
            block = v.read([(kblock-1)*n_per_block*tbin+1 kblock*n_per_block*tbin], 'native');
            block = brick.bin(block, [xbin xbin 1 tbin]);
            a(:,:,:,(kblock-1)*n_per_block+1:kblock*n_per_block) = permute(block, [2 1 3 4]);
        end
    end
% catch
%     if ~isempty(frames), frames = {frames}; end
%     a = aviread(filename,frames{:});
%     switch size(a(1).cdata,3)
%         case 1
%             a = cat(3,a.cdata);
%         case 3
%             s = size(a(1).cdata);
%             nt = length(a);
%             a = cat(2,a.cdata);
%             a = reshape(a,[s(1) s(2) nt 3]);
%         otherwise
%             error('problem')
%     end
% end
if dopermute
    if dodisplay, disp 'transposing frames', end
    a = permute(a,[2 1 3 4]);
end