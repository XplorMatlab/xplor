function xpercent
% function xpercent
%---
% Add % sign to x ticks
%
% See also brick.ypercent, brick.xypercent, brick.ticks

set(gca,'xticklabel',brick.num2str(get(gca,'xtick'),'%i%%','cell'))
