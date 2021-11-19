function V = xplor(data, varargin)
% XPLOR, dynamic multidimensional data viewer
% ---
% function V = xplor(data, options...)
% ---
%
% Input:
% - data        data to visualize: a Matlab ND array or an xplr.XData
%               object
% - options     name-value pairs
%
% Available options:
% - 'header'    header description of the data, either an xplr.Header
%               object of length the number of dimensions of data or a
%               cell array description (see xplr.XData for the syntax)
% - 'name'      name to give to the data
% - 'view'      select which dimensions to view (others will be filtered)
% - 'ROI'       select dimensions in which to select regions of interest
% - 'view&ROI'  select dimensions to view and select ROIs in these
%               dimensions
% - 'display mode'      'time courses' or 'image'
% - 'colormap'          n*3 array or the name of a recognized color map
%                       (this automatically set display mode to 'image')
% - 'controls'          'on'/True (default) or 'off'/False - show/hide the control panel
% ---
% simply type 'xplor' to launch the data import wizard
% type 'xplor demo' to select a range of demos
% type 'xplor test' to launch the "XPLOR logo" demo
%
% See also xplr.XData

% When deployed, make sure the output is defined otherwise an error will
% occur
if isdeployed
    V = [];
end

% Init the log
xplr.log_to_file([], 'init')

% Check license when deployed
if isdeployed && ~optimage_checklicense()
    return
end

% Lauch a demo if no argument
if nargin == 0
    xplr.wizard;
    return
elseif nargin == 1 && ischar(data) 
    switch data
        case 'demo'
            xplor.demo
            return
        case 'test'
            xplor.demo.logo
            return
    end
end

% % drag & drop window
% if nargin==0
%     xplor_window()
%     return
% end

% Gather options inside a structure
for i = 2:2:length(varargin), varargin{i} = {varargin{i}}; end %#ok<CCAT1>
options = reshape(varargin, 2, []);
options(1, :) = brick.strrep(options(1, :), ' ', '_', '&', '_and_');
options = struct(options{:});

% Convert Matlab array to xplr.XData
% create headers and launch view
% if a parameter data is present -> evaluate it (execute and get result)
if isa(data, 'xplr.XData')
    name = data.name;
elseif ischar(data)
    name = data;
    data = evalin('base', name);
    % if no parameters are presents (same case as previous create a demo data)
elseif isfield(options,'name')
    name = options.name;
    options = rmfield(options, 'name');
else
    % no name defined
    name = '';
end
if isnumeric(data) && isempty(data)
    error 'Data is empty'
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

% check options
option_names = fieldnames(options);
valid_names =  {'controls', 'colormap', 'view', 'ROI', 'view_and_ROI', ...
    'filter', 'displaymode', 'display_mode'};
invalid_names = setdiff(option_names, valid_names);
if ~isempty(invalid_names)
    error(['invalid XPLOR option(s): ' brick.strcat(invalid_names, ', ')])
end

% create view
V = xplr.View(data, options);

% apply options
if isfield(options, 'controls')
    V.control_visible = brick.boolean(options.controls);
end
if isfield(options, 'colormap')
    V.D.color_map.c_map_def = options.colormap;
end
actions = {'view', 'ROI', 'view_and_ROI', 'filter'};
for i = 1:length(actions)
    name = actions{i};
    if isfield(options, name)
        action_name = brick.switch_case(name, ...
            'filter', 'add_filter', ...
            name);
        V.C.dim_action(action_name, options.(name))
    end
end
if isfield(options, 'display_mode')
    V.D.display_mode = options.display_mode;
elseif isfield(options, 'displaymode')
    V.D.display_mode = options.displaymode;
end

%---
function xplor_window

hf = brick.figure('[XPLOR]', [220, 120]);
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
