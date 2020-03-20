function V = xplor(data,varargin)
% function V = xplor(data)
% function V = xplor()         [demo data]
% ---
% xplor starter
%
% Input:
% - data        data to visualize 
% - varargin    1-by-N cell array, where N is the number of inputs that
% the function receives after the parameter data

% create a demo xdata with header
if nargin==0
    [~,~,~,dat] = flow;
    dat = permute(dat,[1 3 2]);
    s = size(dat);
    head = xplr.header({'x' 'px' s(1)},{'y' 'px' s(2)},{'axis' num2cell(['a':'y' 'A':'Y'])});
    data = xplr.xdata(dat,head,'Flow Data');
end

% create headers and launch view
% if a parameter data is present -> evaluate it (execute and get result)
if ischar(data)
    name = data;
    data = evalin('base',name);
% if no parameters are presents (same case as previous create a demo data)
elseif nargin==0
    name = 'demo data';
% TODO: document this case
else
    name = inputname(1);
end
if ~isa(data,'xplr.xdata')
    data = squeeze(data); % remove singleton dimensions
    % does user expect an output?
    if nargout > 0
        header = xplr.editHeader(data);
        V = launch_view({data, header, name}, varargin{:});
    else
        xplr.headerEdit(data, @(header)launch_view({data, header, name}, varargin{:}));
        V = [];
    end
else
    if isempty(data.name), data.setName(name); end
    V = launch_view(data, varargin{:});
end
if nargout == 0, clear V, end

%---
function V = launch_view(x, varargin)

% combine data and header
if iscell(x)
    data = xplr.xdata(x{1}, x{2}, x{3});
else
    data = x;
end

% go!
V = xplr.view(data,varargin{:});
