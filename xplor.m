function V = xplor(data, varargin)
% function V = xplor(data, options...)
% function V = xplor()         [demo data]
% ---
% xplor starter
%
% Input:
% - data        data to visualize: a Matlab ND array or an xplr.XData
%               object
% - options     name-value pairs
%
% Available options:
% - 'header'    header description of the data, either an xplr.Header
%               object of length the number of dimensions of data or a
%               description (see xplr.XData for the syntax)
% - 'name'      name to give to the data
% - 'view'      select which dimensions to view (others will be filtered)

% Create a demo xdata with header
if nargin == 0
    [~, ~, ~, dat] = flow;
    dat = permute(dat, [1, 3, 2]);
    s = size(dat);
    head = xplr.Header( ...
        {'x', 'px', s(1)}, {'y', 'px', s(2)}, ...
        {'axis', num2cell(['a':'y' 'A':'Y'])} ...
        );
    data = xplr.XData(dat, head, 'Flow Data');
end

% % drag & drop window
% if nargin==0
%     xplor_window()
%     return
% end

% Gather options inside a structure
for i = 2:2:length(varargin), varargin{i} = {varargin{i}}; end %#ok<CCAT1>
options = struct(varargin{:});

% Convert Matlab array to xplr.XData
% create headers and launch view
% if a parameter data is present -> evaluate it (execute and get result)
if ischar(data)
    name = data;
    data = evalin('base', name);
    % if no parameters are presents (same case as previous create a demo data)
elseif isfield(options,'name')
    name = options.name;
    options = rmfield(options, 'name');
else
    name = inputname(1);
end
if isfield(options,'header')
    if isa(data, 'xplr.XData')
        error 'header description in the options will be ignored since data is already an xplr.XData object'
    end
    data = xplr.XData(data,options.header, name);
    options = rmfield(options, 'header');
end
if isa(data, 'xplr.XData')
    if isempty(data.name), data.set_name(name); end
    V = launch_view(data, options);
else
    data = squeeze(data); % remove singleton dimensions
    % does user expect an output?
    if nargout > 0
        % wait for xplr.edith_header to finish
        header = xplr.edit_header(data);
        V = launch_view({data, header, name}, options);
    else
        % do not wait, view will be launched later
        xplr.HeaderEdit( ...
            data, ...
            @(header)launch_view({data, header, name}, ...
            options) ...
            );
        V = [];
    end
end
if nargout == 0, clear V, end

%---
function V = launch_view(x, options)

% combine data and header
if iscell(x)
    data = xplr.XData(x{1}, x{2}, x{3});
else
    data = x;
end

% create view
V = xplr.View(data);

% apply options
option_names = fieldnames(options);
for i = 1:length(option_names)
    name = option_names{i};
    value = options.(name);
    switch name
        case 'view'
            V.C.dim_action('viewdim', value)
        otherwise
            disp(['invalide xplor option ''' name ''''])
    end
end

%---
function xplor_window

hf = fn_figure('[XPLOR]', [220, 120]);
delete(findall(hf, 'parent', hf))

uicontrol('unit', 'normalized', 'pos', [0 .6 1 .4], ...
    'style', 'pushbutton', 'enable', 'inactive', 'fontsize', 20, ...
    'backgroundcolor', [.8, .8, 1], ...
    'string', 'XPLOR', ...
    'callback', 'disp hello');
uicontrol('unit', 'normalized', 'pos', [0, 0, 1, .6], ...
    'style', 'edit', 'max', 2, 'fontsize', 14, ...
    'string', 'Drag & Drop variable here to XPLOR it', ...
    'ForegroundColor', [1, 1, 1]*.5, 'fontangle', 'italic', ...
    'callback', 'disp hello');

%     [~,~,~,dat] = flow;
%     dat = permute(dat,[1 3 2]);
%     s = size(dat);
%     head = xplr.header({'x' 'px' s(1)},{'y' 'px' s(2)}, ...
%       {'axis' num2cell(['a':'y' 'A':'Y'])});
%     data = xplr.xdata(dat,head,'Flow Data');
