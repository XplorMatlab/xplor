function [x, clip] = clip(x,varargin)
%brick.clip Rescale data, restrict the range, color
%---
% function [x clip] = clip(x[,clipflag][,outflag][,nanvalue])
%---
% Rescale and restrict to a specific range ("clip") the data in an array.
% Make a color image if requested.  
%
% Input:
% - x           array (any dimension)
% - clipflag    clipping mode:
%               [a b]                   define manually min and max value
%               'fit','mM' or 'minmax'  use minimum and maximum [default]
%               'Xstd'                  use mean and X times standard deviation
%               'prcA-B'                use percentiles (if B is omitted,
%                                       use B = 100-A; if B<30, uses 100-B)
%               add '[value]' at the end (e.g. 'fit[0]') to center the
%               clipping range on the specified value 
%               add '[value|position]' at the end (e.g. 'fit[0|.1]') to fix
%               a given value (here 0) at a given position (between 0 and
%               1, here .1 -> 0 value will be kept at bottom)
%               
% - outflag     output format
%               [a b]       define minimum and maximum value [default, with
%                           a=0 and b=1]
%               n           integer values between 1 and n (output will be
%                           of class uint8 if n<=256, uint16 if n<=65536,
%                           etc.
%               'uint8', 'uint16', ..   integer values between 0 and max
%               nx3 array   returns a (n+1)-dimensional array using this
%                           colormap 
%               char array  use this colormap (for example 'jet' -> use
%                           jet(256))
%               'scaleonly' rescale data but do not coerce within the range
%               'getrange'  output not the cliped data, but the calculated
%                           clipping range (can be useful for 'std' and
%                           'prc' clipping calculations)
%
% - nanvalue    value to give to NaNs
%
% Output:
% - x           the clipped image (or the clipping range if outflag is
%               'getrange') 
% - clip        the clipping range 

% Thomas Deneux
% Copyright 2007-2017

if nargin==0, help brick.clip, return, end

% Input
xc = x(:);  % column vector
clipflag = []; outflag = []; nanvalue = [];
for k=1:length(varargin)
    a = varargin{k};
    if ischar(a)
        if any(regexp(a,'^fit|mM|minmax')) || any(regexpi(a,'(^prc)|((st|sd|std))'))
            clipflag = a;
        elseif regexp(a,'^[0-9e\-.]+ +[0-9e\-.]+$') % two numbers
            clipflag = str2num(a); %#ok<ST2NM>
        else
            outflag = a;
        end
    else
        if isvector(a) && length(a)==2 && isempty(clipflag)
            clipflag = a;
        elseif isempty(outflag)
            outflag = a;
        else
            nanvalue = a;
        end
    end
end
if isempty(clipflag), clipflag='mM'; end
if isempty(outflag), outflag=[0 1]; end

% Compute clipping range
if isnumeric(clipflag)
    if ~isvector(clipflag) || length(clipflag)~=2, error('clipping vector must have 2 elements'), end
    clip = brick.row(clipflag);
else
    ibaseline = regexp(clipflag,'\[.*\]$');
    if isempty(ibaseline)
        base_value=[];
    else
        [base_value, base_position] = brick.regexptokens(clipflag(ibaseline+1:end-1), '([^\|]*)\|?([^\|]*)');
        base_value = str2double(base_value);
        if isempty(base_position)
            base_position = .5;
        else
            base_position = str2double(base_position);
        end
        clipflag = clipflag(1:ibaseline-1);
    end
    xstd = regexpi(clipflag,'^([\d.]*)(st|sd|std)$','tokens');
    if ~isempty(xstd)
        xstd = xstd{1}{1};
    else
        xstd = regexpi(clipflag,'^(st|sd|std)([\d.]*)$','tokens');
        if ~isempty(xstd), xstd = xstd{1}{2}; end
    end
    xprc = regexpi(clipflag,'^prc([\d.]*)[-_]*([\d.]*)$','tokens');
    if brick.ismemberstr(clipflag,{'fit' 'mM' 'minmax'})
        if isempty(base_value)
            clip = [min(xc) max(xc)];
        else
            M = max(xc) - base_value;
            m = base_value - min(xc);
            if base_position == 0
                % ignore values below baseline
                clip = base_value + [0 M];
            elseif base_position == 1
                % ignore values above baseline
                clip = base_value + [-m 0];
            else
                e = max(M/(1-base_position), m/base_position);
                clip = base_value + [-base_position 1-base_position]*e;
            end
        end
    elseif ~isempty(xstd)
        xc = brick.float(xc);
        if isempty(xstd), xstd=1; else xstd=str2double(xstd); end
        if isempty(base_value)
            m = brick.nmean(xc); 
            st = brick.nstd(xc);
            base_position = .5;
        else
            % use baseline value instead of mean to compute "standard deviation"
            m = base_value; 
            st = sqrt(brick.nmean((xc-m).^2));
        end
        clip = m + 2*[-base_position 1-base_position]*xstd*st;
    elseif ~isempty(xprc)
        low = str2double(xprc{1}{1});
        high = str2double(xprc{1}{2});
        if isempty(base_value)
            if isnan(high), high=100-low; elseif high<30, high=100-high; end
            clip = [brick.prctil(xc,low) brick.prctil(xc,high)];
        else
            if base_position == 0
                % ignore values smaller than baseline
                if isnan(high), high=100-low; elseif high<30, high=100-high; end
                clip = [base_value brick.prctil(xc(xc>=base_value),high)];
            elseif base_position == 1
                % ignore values higher than baseline
                clip = [brick.prctil(xc(xc<=base_value),low) base_value];
            else
                if ~isnan(high)
                    warning 'cannot set independently the percentile of low and high out-of-range when center value is fixed'
                    if high<30, high=100-high; end
                    low = (low+high)/2;
                end
                above = (xc >= base_value);
                dev = xc;
                dev(above) = (xc(above)-base_value) / (1-base_position);
                dev(~above) = (base_value-xc(~above)) / base_position;
                max_dev = brick.prctil(dev, 100-low);
                clip = base_value + [-base_position 1-base_position]*max_dev;
            end
        end
    else
        error('erroneous clipping option')
    end
end
clip = full(clip); % if x was sparse, clip was sparse so far
if diff(clip)==0, clip = double(clip)+[-1 1]; end

% Check output mode
doclip = true;
if strcmp(outflag,'getrange')
    x = clip;
    return
elseif strcmp(outflag,'scaleonly')
    doclip = false;
elseif ischar(outflag) && any(strfind(outflag,'uint'))
    docolor = false;
    n = intmax(outflag);
elseif ischar(outflag)
    docolor = true;
    fname = outflag;
    if strfind(fname,'.LUT')
        cm = brick.readasciimatrix(fname);
    else
        cm = feval(fname,256);
    end
    [n, nc] = size(cm);
elseif isscalar(outflag)
    docolor = false;
    n = outflag;
    if mod(n,1) || n<=0, error('scalar for output format must be a positive integer'), end
elseif isvector(outflag) && length(outflag)==2
    docolor = false;
    n = 0;
    a = outflag(1);
    b = outflag(2);
else
    nc = size(outflag,2);
    docolor = true;
    cm = outflag;
    n = size(cm,1);
end

% Convert data to float if the data is integer (we avoided doing it before
% because it was not needed in the case of the 'getrange' output flag)
[x, clip] = deal(brick.float(x),brick.float(clip));

% Clip data
x = (x-clip(1))/diff(clip);
if ~doclip, return, end
if ~isempty(nanvalue)
    xnan = isnan(x);
end
if isa(x,'double')
    upperbound = 1-eps(1); % it is convenient that 1 cannot be reached
else
    upperbound = 1-eps(single(1));
end
x = min(upperbound,max(0,x)); 

% Scaling
if n
    if isinteger(n)
        x = cast(floor(double(n)*x),'like',n); % values between 0 and n-1
    else
        x = floor(n*x)+1; % values between 1 and n
    end
elseif a==0 && b==1
    % nothing to do
else
    x = a + x*(b-a); % b cannot be reached
end

% Color and value for NaNs
if docolor
    s = size(x); %if s(2)==1, s(2)=[]; end
    if ~isempty(nanvalue)
        x(xnan) = n+1;
        cm(n+1,:) = nanvalue;
    end
    x = reshape(cm(x,:),[s nc]);
    if length(s)>2, x = permute(x,[1 2 length(s)+1 3:length(s)]); end
else
    if ~isempty(nanvalue)
        x(xnan) = nanvalue;
    end
end




        
