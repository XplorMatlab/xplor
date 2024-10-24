function markpvalue(x,y,p,varargin)
%MARKPVALUE Draw stars to mark significancy of results
%---
% function markpvalue(x,y,p[,'ns|p'][,'vertical'][,text options...])
%---
% Mark p-value if significant with star(s) at location x,y. 
% If x and y are 2-elements vectors, draws a line segment (x,y(1)) and mark
% star(s) at (mean(x),y(2)). If y is a 3-element vector, adds also small
% vertical limits.

% Thomas Deneux
% Copyright 2015-2017

if isnan(p), return, end

% Input
topt = varargin;
kparent = find(strcmp(varargin(1:2:end),'parent'));
if ~isempty(kparent)
    popt=varargin(2*kparent-1:2*kparent); ha = popt{2};
    topt(2*kparent-1:2*kparent)=[];
else 
    popt={}; ha = gca;
end
displaymode = 'default'; orientation = 'horizontal';
while ~isempty(topt) && ischar(topt{1})
    switch lower(topt{1})
        case {'ns' 'all'}
            displaymode = lower(topt{1});
        case 'p'
            displaymode = 'pvalue';
        case 'vertical'
            orientation = 'vertical';
        otherwise
            break
    end
    topt(1) = [];
end

% build string to display
fontweight = 'normal';
if ischar(p)
    stars = p;
    if contains(p, '<')
        fontweight = 'bold'; 
    end
elseif strcmp(displaymode,'pvalue')
    stars = num2str(p,'p=%.2g');
    if p<=.05, fontweight = 'bold'; end
else
    doNS = ismember(displaymode,{'ns' 'all'});
    nstar = floor(log10(1/p));
    if p>.05
        if doNS, stars = 'n.s.'; else return, end
    elseif p==0
        stars = '*!';
    elseif nstar<=5
        stars = brick.switch_case(orientation, ...
            'horizontal',repmat('*',1,nstar), ...
            'vertical',repmat({'*'},1,nstar));
    else
        stars = brick.switch_case(orientation, ...
            'horizontal',['*' num2str(nstar)], ...
            'vertical',{'*' num2str(nstar)});
        %stars = ['p<1e-' num2str(nstar)];
    end
    if strcmp(displaymode,'all')
        stars = [stars num2str(p,' (p=%.3g)')];
    end
end

% display
if isempty(y), ylim = get(ha,'ylim'); y = ylim(1)*.1+ylim(2)*.9; end
xs = mean(x); ys = y(end);
h = text(xs,double(ys),stars,'fontweight',fontweight, ...
    'horizontalalignment','center','verticalalignment','middle',popt{:});
if length(x)==2
    if length(y)==2
        h(2) = line(x,y(1)*[1 1],'color','k',popt{:});
    elseif length(y)==3
        h(2) = line(x([1 1 2 2]),y([1 2 2 1]),'color','k',popt{:});
    else
        error 'incompatible lengths for x and y'
    end
end
brick.set(h,topt{:})
	
