function ypercent
% function ypercent
%---
% Add % sign to y ticks
%
% See also brick.xpercent, brick.xypercent, brick.ticks

set(gca,'yticklabel',brick.num2str(get(gca,'ytick'),'%i%%','cell'))
