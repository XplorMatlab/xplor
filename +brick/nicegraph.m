function nicegraph(ha)
%NICEGRAPH Improve aspect of graph display
%---
% function nicegraph([ha])
%---
% see also brick.labels, brick.scale, brick.plotscale

% Thomas Deneux
% Copyright 2005-2017

if nargin==0, ha = gca; end

ha = findobj(ha,'type','axes');
for hak = brick.row(ha)
    set(hak,'tickdir','out','ticklength',[.03 1],'box','off')
end
