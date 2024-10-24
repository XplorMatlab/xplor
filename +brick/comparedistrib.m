function [p, hl] = comparedistrib(x, y, varargin)
%COMPAREDISTRIB Perform a nonparametric test and display data points and results
%---
% function [pval hl] = comparedistrib(x,y[,test][,'tail','left|right|both']
%       [,'showmean'][,'ylim',ylim][,'xlabels',xlabels][,'pdisplaymode','ns|p']
%       [,'jitter|xjittrand|xdispatch'[,jittersize]],['xpos',[x1, x2]]
%       [,'marker',marker][,'markersize',markersize])
%---
% Perform any of 'ranksum', 'signrank' or 'signtest' test and display the
% data and p-value.
%
% Input
% - x,y     data points; for signrank or signtest, y can be a scalar
%           (the tested median/mean value, typically 0)
% - test    'ranksum' (=default if y is nonscalar)
%           'signrank'
%           'signtest' (=default if y is scalar)
%           'bootstrap' (test on the mean)
%           'bootstrapmedian' (test on the median)
%           'bootstrapsign', 'bootstrapsignmedian'
%           p - providing a p-value results in skipping the test and
%           displaying this p-value
% - 'jitter', 'xjittrand', 'xdispatch'  different modes for slightly moving
%           points either vertically ('jitter') or horizontally, randomly
%           or not ('xdispatch'), so as to see their multiplicity;
%           'xdispatch' will make appear a distribution vertically
%
% See also brick.markpvalue

% Thomas Deneux
% Copyright 2015-2017

% Input: data
if nargin<2
    if isvector(x)
        y = 0;
    else
        if size(x,2)~=2
            error 'single matrix input must have two columns'
        end
        [x, y] = deal(x(:,1),x(:,2));
    end
end

% Input: options
i = 0; tail = 'both'; ylim = []; showmean = false; xlabels = {}; method = [];
xpos = [1 2]; marker = 'o'; markersize = 'default';
pdisplaymode = 'ns'; jittermode = ''; jittersize = [];
while i<length(varargin)
    i = i+1;
    switch(varargin{i})
        case 'tail'
            i = i+1;
            tail = varargin{i};
        case {'left', 'right'}
            tail = varargin{i};
        case 'xlabels'
            i = i+1;
            xlabels = varargin{i};
        case 'ylim'
            i = i+1;
            ylim = varargin{i};
        case 'pdisplaymode'
            i = i+1;
            pdisplaymode = varargin{i};
        case 'showmean'
            showmean = true;
        case {'signtest' 'ranksum' 'signrank' ...
                'bootstrap' 'bootstrapmedian' 'bootstrapsign' 'bootstrapsignmedian' ...
                'anova'}
            method = varargin{i};
        case {'p' 'ns'}
            pdisplaymode = varargin{i};
        case {'jitter', 'xjittrand' 'xdispatch'}
            jittermode = varargin{i};
            if i<length(varargin) && isnumeric(varargin{i+1})
                i = i+1;
                jittersize = varargin{i};
            end
        case 'xpos'
            i = i+1;
            xpos = varargin{i};
        case 'marker'
            i = i+1;
            marker = varargin{i};
        case 'markersize'
            i = i+1;
            markersize = varargin{i};
        otherwise
            error('unknown flag ''%s''',varargin{i})
    end
end
if isempty(method)
    method = brick.switch_case(isscalar(y),'signtest','ranksum');
end

% is population the same in x and y?
same_population = any(strfind(method, 'sign')) && ~isscalar(y);
if same_population
    assert(length(x) == length(y))
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok);
    y = y(ok);
else
    x = x(~isnan(x));
    y = y(~isnan(y));
end
nx = length(x);
ny = length(y);


% p-value
if all(isnan(x)) || all(isnan(y))
    p = NaN;
elseif isnumeric(method)
    p = method;
else
    switch method
        case 'anova'
            groups = [zeros(1, nx) ones(1, ny)];
            p = anova1([brick.row(x) brick.row(y)], groups, 'off');
        case {'ranksum' 'signrank' 'signtest'}
            p = feval(method,x,y,'tail',tail);
        case 'bootstrap'
            p = brick.bootstrap(x,y,'mean','tail',tail,'npermmax',2e5);
        case 'bootstrapmedian'
            p = brick.bootstrap(x,y,'median','tail',tail,'npermmax',2e5);
        case 'bootstrapsign'
            p = brick.bootstrap(x-y,[],'mean','tail',tail,'npermmax',2e5);
        case 'bootstrapsignmedian'
            p = brick.bootstrap(x-y,[],'median','tail',tail,'npermmax',2e5);
        otherwise
            error('unknown test ''%s''',method);
    end
end


% move points in x or y
xjit = {zeros(1,nx), zeros(1,ny)};
switch jittermode
    case ''
        % no jittering
    case 'jitter'
        % add a jitter in y to displayed data to help see the multiplicity
        % of identical points
        if isempty(jittersize)            
            m = min(min(x),min(y));
            M = max(max(x),max(y));
            jittersize = (M-m) / 100;
        end
        x = x + randn(size(x)) * jittersize;
        y = y + randn(size(y)) * jittersize;
    case 'xjittrand'
        % move points randomly in x so that they will reflect the
        % distribution, which will appear as a vertical curve
        if isempty(jittersize), jittersize = .3; end
        X = {x, y};
        for k = 1:2
            xk = X{k};
            [counts, edges] = histcounts(xk, 100);
            if all(counts <= 1)
                continue
            end
            curve = (counts-1) / max(counts-1);
            curve = interp1( ...
                [edges(1) (edges(1:end-1)+edges(2:end))/2 edges(end)], ...
                curve([1 1:end end]), ...
                xk);
            xjit{k} = (2*rand(1, length(xk),1)-1) .* curve * jittersize;
        end
    case 'xdispatch'
        % move points in x in a deterministic manner so that we will
        % see each individual point, and the distribution will appear as a
        % vertical curve
        X = {x, y};
        for k = 1:2
            xk = X{k};
            range = max(xk) - min(xk);
            if range == 0
                c = ones(1, length(xk));
                nval = 1;
            else
                nbin = 35;
                [u, ~, c] = unique(round(xk / (max(xk)-min(xk)) * nbin));
                c = c(:)'; % row vector
                nval = length(u);
            end
            if nval == length(xk)
                % all values are unique
                continue
            end
            counts = zeros(1, nval);
            pos = zeros(1, nval);
            for i = 1:length(c)
                ci = c(i);
                pos(i) = counts(ci);
                counts(ci) = counts(ci)+1;
            end
            xjit{k} = pos - (counts(c)-1)/2;
        end
        if isempty(jittersize)
            M = max(max(xjit{1}), max(xjit{2}));
            if M == 0
                jittersize = 0;
            else
                jittersize = .3 / M;
            end
        end
        for k = 1:2
            xjit{k} = xjit{k} * jittersize;
        end
    otherwise        
        error('invalid points mode ''%s''', jittermode)
end

% display
singledistrib = isscalar(y) && ~strcmp(method,'ranksum');
if singledistrib
    xpos = xpos(1);
end
alldata = [brick.row(x) brick.row(y)];
if strcmp(method,'ranksum')
    % no connecting lines 
    % group points according to whether they are below/above the median
    [~, ord] = sort(x);
    idx = ord(1:floor(nx/2));
    a = plot(xpos(1)*ones(1,floor(nx/2))+xjit{1}(idx), ...
        x(idx),marker,'color',[1 1 1]*.5,'markersize',markersize);
    hold on
    idx = ord(floor(nx/2)+1:nx);
    b = plot(xpos(1)*ones(1,ceil(nx/2))+xjit{1}(idx), ...
        x(idx),marker,'color',[1 1 1]*.6,'markersize',markersize);
    [~, ord] = sort(y);
    idx = ord(1:floor(ny/2));
    c = plot(xpos(2)*ones(1,floor(ny/2))+xjit{2}(idx), ...
        y(idx),marker,'color',[1 1 1]*.5,'markersize',markersize);
    idx = ord(floor(ny/2)+1:ny);
    d = plot(xpos(2)*ones(1,ceil(ny/2))+xjit{2}(idx), ...
        y(idx),marker,'color',[1 1 1]*.6,'markersize',markersize);
    hold off
    hl{1} = [a b c d];
elseif ~same_population
    % no connecting lines 
    a = plot(xpos(1)*ones(1,nx)+xjit{1}, ...
        x,marker,'color',[1 1 1]*.6,'markersize',markersize);
    if singledistrib
        hl{1} = a;
    else
        hold on
        b = plot(xpos(2)*ones(1,ny)+xjit{2}, ...
            y,marker,'color',[1 1 1]*.6,'markersize',markersize);
        hold off
        hl{1} = [a b];
    end
else
    hl{1} = plot([xpos(1)+xjit{1}; xpos(2)+xjit{2}], ...
        [brick.row(x); brick.row(y)], ...
        'color',[1 1 1]*.6, ...
        'marker',marker,'markersize',markersize); % connecting lines
end

% display means and/or medians
if showmean
    if singledistrib
        line(xpos + [-.3 .3], brick.nmean(x) * [1 1], 'color','b')
    else
        line(xpos,[brick.nmean(x) brick.nmean(y)],'color','b')
    end
end
switch method
    case {'ranksum' 'bootstrapmedian'}
        % show individual medians
        if singledistrib
            hl{2} = line(xpos,brick.nmedian(x),'color','k','linestyle','none','marker','*');
        else
            hl{2}(1) = line(xpos,[brick.nmedian(x) brick.nmedian(y)],'color','k','linestyle','none','marker','*');
            hl{2}(2) = line(xpos,[brick.nmedian(x) brick.nmedian(y)],'color','k','linewidth',2);
        end
    case {'anova' 'bootstrap'}
        % show individual means
        if singledistrib
            hl{2} = line(xpos,brick.nmean(x),'color','k','linestyle','none','marker','*');
        else
            hl{2}(1) = line(xpos,[brick.nmean(x) brick.nmean(y)],'color','k','linestyle','none','marker','*');
            hl{2}(2) = line(xpos,[brick.nmean(x) brick.nmean(y)],'color','k','linewidth',2);
        end
    otherwise
        % show individual means (not medians), but also a slope indicating the
        % median difference (which is different from the difference
        % of the medians!)
        if singledistrib
            hl{2} = line(xpos,brick.nmean(x),'color','k','linestyle','none','marker','*');
        else
            hl{2}(1) = line(xpos,[brick.nmean(x) brick.nmean(y)],'color','k','linestyle','none','marker','*');
            yl = mean([brick.nmedian(x) brick.nmedian(y)])+[-.5 .5]*brick.nmedian(y-x);
            hl{2}(2) = line(xpos,yl,'color','k','linewidth',2);
        end
end

% limits
if singledistrib
    xlim = xpos + [-1 1];
else
    xlim = xpos + [-1 1]*diff(xpos);
end
if isempty(ylim)
    if all(isnan(alldata))
        ylim = [0 1];
    else
        m = min(alldata); M = max(alldata);
        ylim = m+[-.1 1.3]*(M-m);
    end
end
set(gca,'xlim',xlim,'ylim',ylim)

% p-value
if contains(method, 'bootstrap') && p < 1e-5
    p = 'p<1e-5';
end
brick.markpvalue(mean(xpos),[],p,pdisplaymode)

% x ticks
if ~isempty(xlabels)
    set(gca,'xtick',1:length(xlabels),'xticklabel',xlabels,'xTickLabelRotation',30)
end

% reference
if singledistrib
    h = line(xlim, [y y], 'linestyle', '--', 'color', 'k');
    uistack(h, 'bottom')
end

% output?
if nargout==0
    clear p
end

% immediate display is usefull when multiple comparisons are being computed
drawnow
