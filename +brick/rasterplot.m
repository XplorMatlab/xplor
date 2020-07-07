function [hl hlsep] = rasterplot(times,varargin)
%RASTERPLOT Display multiple sets of events as small bars
%---
% function [hlspikes hlsep] = rasterplot(times,[other arguments as for spikedisplay])
%---
% This is a wrapper of brick.spikedisplay. Additionally to calling
% brick.spikedisplay, it initially clears the axes and set a new axis.
%
% See also brick.spikedisplay

% Thomas Deneux
% Copyright 2008-2017

if nargin==0, help brick.rasterplot, help brick.spikedisplay, return, end

iparent = find(strcmp(varargin,'parent'));
if ~isempty(iparent)
    ha = varargin{iparent+1};
else
    ha = gca;
end

% clear axes
cla(ha)

% auto-pos?
if ~iscell(times), times = {times}; end
if nargin<2 || ischar(varargin{1})
    if ~isvector(times)
        [n1 n2] = size(times); 
        y = brick.add(-.4+(.5:n1)'/n1*.8,1:n2);
        spikeheight = .8/n1*.7;
    else
        n = length(times);
        y = 1:n;
        spikeheight = .6;
    end
    varargin = [{y spikeheight} varargin];
end

% display spikes
[hl y] = brick.spikedisplay(times,varargin{:});

% decorate
axis(ha,'tight')
ax = axis;

% separate traces
if ~isvector(times)
    box(ha,'on')
    ngroup = size(times,2);
    y = mean(y,1); % one value per column
    hlsep = zeros(1,ngroup-1);
    for i=1:ngroup-1
        hlsep(i) = line(ax(1:2),[1 1]*mean(y([i i+1])),'color','k');
    end
    set(gca,'ytick',1:ngroup)
end

% output
if nargout==0, clear hl, end
