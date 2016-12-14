function V = xplor(data,varargin)

% create a demo data
if nargin==0
    [~,~,~,dat] = flow;
    dat = permute(dat,[1 3 2]);
    s = size(dat);
    head = xplr.header({'x' 'px' s(1)},{'y' 'px' s(2)},{'axis' num2cell(['a':'y' 'A':'Y'])});
    data = xplr.xdata(dat,head,'Flow Data');
end

% create headers
if ~isa(data,'xplr.xdata')
    if ischar(data)
        name = data;
        data = evalin('base',name);
    else
        name = inputname(1);
    end
    data = squeeze(data); % remove singleton dimensions
    head = xplr.editHeader(data);
    if isempty(head)
        % header editing was canceled
        if nargout, V = []; end
        return
    end
    data = xplr.xdata(data,head,name);
end

% go!
V = xplr.view(data,varargin{:});
if nargout==0, clear V, end