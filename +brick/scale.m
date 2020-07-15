function hx = scale(barsize,label,col)
%SCALE Scale bar for image display
%---
% function hx = scale(barsize,label[,color])
%---
% Remove ticks from image display and draw a scale bar with the size and
% label as specified
% 
% See also brick.plotscale, brick.nicegraph, brick.labels

% Thomas Deneux
% Copyright 2007-2017

if nargin==0, help brick.scale, return, end
if nargin<3, col='black'; end    

set(gca,'xtick',[],'ytick',[])
delete(findobj(gca,'tag','brick.scale'))

if ~strcmp(get(gca,'DataAspectRatioMode'),'manual'), axis image, end
ax = axis; 

% bar starts at (10,10) right and above bottom-left corner and has length
% 'barize'; the text is 10 pixels above
lineorigin = brick.coordinates('b2a',[15 15]','position');
xdir = 1-2*strcmp(get(gca,'xdir'),'reverse');
lineend = lineorigin + xdir*[barsize 0]';
linepos = [lineorigin lineend];
textpos = mean(linepos,2) + ...
    brick.coordinates('b2a',[0 10]','vector');

hx(1) = line(linepos(1,:),linepos(2,:),'color',col,'linewidth',3, ...
    'tag','brick.scale');
hx(2) = text(textpos(1),textpos(2),label, ...
    'horizontalalignment','center','verticalalignment','middle', ...
    'color',col,'tag','brick.scale');
if nargout==0, clear hx, end
