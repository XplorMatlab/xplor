function V = xplor(data,varargin)
% function V = xplor(data,varargin)
    % function V = xplor(data)
    % function V = xplor()         [demo data]
    %---
    % Base xplor starter
    %
    % Input:
    % - data        data to visualize 
    % - varargin    number of arguments in input
    %
    % 
% create a demo data
if nargin==0
    [~,~,~,dat] = flow; 
    dat = permute(dat,[1 3 2]);
    s = size(dat);
    head = xplr.header({'x' 'px' s(1)},{'y' 'px' s(2)},{'axis' num2cell(['a':'y' 'A':'Y'])});
    data = xplr.xdata(dat,head,'Flow Data');
end

% create headers and launch view
if ischar(data)
    name = data;
    data = evalin('base',name);
elseif nargin==0
    name = 'demo data';
else
    name = inputname(1);
end
if ~isa(data,'xplr.xdata')
    data = squeeze(data); % remove singleton dimensions
    xplr.headerEdit(data, @(header)launch_view({data, header, name}, varargin{:}));
else
    launch_view(data, varargin{:})
end

%---
function launch_view(x, varargin)

% combine data and header
if iscell(x)
    data = xplr.xdata(x{1}, x{2}, x{3});
else
    data = x;
end

% go!
xplr.view(data,varargin{:});
