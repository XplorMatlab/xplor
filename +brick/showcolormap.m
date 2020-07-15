function ha = showcolormap(varargin)
%SHOWCOLORMAP Display a color map in a given axes or in a separate figure
%---
% function [ha =] showcolormap([cmap][,clip][,ha|'autopos'])
%---
% opens a new figure and display in it the colormap at the given clipping
% value

% Thomas Deneux
% Copyright 2011-2017


% input
cmap = []; clip = []; ha = []; doautopos = false;
for k=1:length(varargin)
    a = varargin{k};
    if ischar(a) && strcmp(a,'autopos')
        doautopos = true;
    elseif ischar(a) || size(a,2)==3
        cmap = a;
        if ischar(cmap), cmap = feval(cmap,256); end
    elseif isscalar(a) && ishandle(a)
        ha = a;
    else
        clip = a;
    end
end
if isempty(cmap), cmap = get(gcf,'colormap'); end
if isempty(clip)
    if isempty(findobj(0,'type','axes'))
        clip = [0 1];
    else
        clip = get(gca,'clim');
    end
end
if doautopos
    ha0 = gca;
    try delete(getappdata(ha0,'brick_showcolormap')), end %#ok<TRYNC>
    hf = get(ha0,'parent');
    units = get(ha0,'units');
    set(ha0,'units','normalized');
    p = get(ha0,'pos');
    figsiz = brick.pixelsize(hf);
    w = min((1-(p(1)+p(3)))/5,22/figsiz(1)); h = p(4)*2/3;
    ha = axes('pos',[p(1)+p(3)+3*w p(2)+h/4 w h]);
    set([ha0 ha],'units',units)
    setappdata(ha0,'brick_showcolormap',ha)
    verticalaxes = true;
elseif isempty(ha)
    hf = figure(873);
    clf(hf)
    brick.setfigsize(hf,200,500);
    set(hf,'color','w','tag','color bar')
    ha = axes('pos',[.5 .2 .1 .7]);
    verticalaxes = true;
else
    pixsiz = brick.pixelsize(ha);
    verticalaxes = pixsiz(2)>pixsiz(1);
end

% display
if verticalaxes
    im = repmat(permute(cmap,[1 3 2]),[1 10]);
    imagesc([0 1],clip,im,'parent',ha)
    set(ha,'ydir','normal','xtick',[])
else
    im = repmat(permute(cmap,[3 1 2]),[10 1]);
    imagesc(clip,[0 1],im,'parent',ha)
    set(ha,'xdir','normal','ytick',[])
end
%if ~defclip, set(ha,'xtick',[],'ytick',[]), end
axis(ha,'normal')
if doautopos, set(ha,'handlevisibility','off'), end % avoid next call to brick.showcolormap with 'autopos' option to be applied to ha

% output?
if nargout==0, clear ha, end

