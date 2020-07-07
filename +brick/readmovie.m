function a = readmovie(filename,varargin)
% function a = readmovie(filename,frames[,'nodisplay'][,'nopermute'])
%---
% read an avi file and stores it into a 2D+time array 
% (2D+time+channel if color movie)
%
% See also brick.savemovie

% Thomas Deneux
% Copyright 2004-2017

if nargin<1
    filename = brick.getfile('*.avi');
end
if ~exist(filename,'file')
    error('file ''%s'' does not exist',filename)
end
frames = {};
dodisplay = true; dopermute = true;
for k=1:length(varargin)
    a = varargin{k};
    if isnumeric(a)
        frames = a;
    elseif ischar(a)
        switch a
            case 'nodisplay'
                dodisplay = false;
            case 'nopermute'
                dopermute = false;
            otherwise
                error argument
        end
    else
        error argument
    end
end

if dodisplay, disp 'reading', end
try
    % recent Matlab version
    if ~isempty(frames), frames = {frames([1 end])}; end
    a = read(VideoReader(filename),frames{:},'native');
catch
    if ~isempty(frames), frames = {frames}; end
    a = aviread(filename,frames{:});
    switch size(a(1).cdata,3)
        case 1
            a = cat(3,a.cdata);
        case 3
            s = size(a(1).cdata);
            nt = length(a);
            a = cat(2,a.cdata);
            a = reshape(a,[s(1) s(2) nt 3]);
        otherwise
            error('problem')
    end
end
if dopermute
    if dodisplay, disp 'transposing frames', end
    a = permute(a,[2 1 3 4]);
end