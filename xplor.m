function V = xplor(varargin)
% XPLOR, dynamic multidimensional data viewer
% ---
% function V = xplor(data, options...)
% function V = xplor(data1, data2, ..., options...)
% ---
%
% Input:
% - data        data to visualize: a Matlab ND array or an xplr.XData
%               object
% - data1, data2, ...   several data will be concatenated along a new
%               dimensions if they are all of the same size, along the
%               dimension where their sizes differe if there is one
% - options     name-value pairs
%
% Available options:
% - 'header'    header description of the data, either an xplr.Header
%               object of length the number of dimensions of data or a
%               cell array description (see xplr.XData for the syntax)
% - 'name'      name to give to the data
% - 'view'      select which dimensions to view (others will be filtered)
%               use a cell array of arrays or cell arrays to specify
%               organization as for 'organize' option below; for example:
%               - {'x' 'y' 'condition' 'day'} or [1 2 3 4] will let xplor
%               determine automatically how to organize dimensions
%               - but {{'x' 'condition'} 'y' 'day'} or {[1 3] 2 4} will
%               specify both dimensions to view and their organization
% - 'filter'    select which dimensions to filter (others will be viewed)
% - 'organize'  specify where the dimensions will appear; should be a 1 to
%               4-elements array or cell array
%               for example:
%               [1 2 3 4]   view dimensions 1 in x-axis, 2 in y-axis, 3 in
%                           a x/y grid arrangement, 4 as superimposed data
%               {'time' 'wavelength'}   view dimension 'time' in x-axis and
%                           'wavelength' in y-axis
%               {{'x' 'condition'} 'y' 'day'}   view dimensions 'x' and
%                           'condition' in x-axis, 'y' in y-axis and 'day'
%                           in x/y grid arrangement
% - 'ROI'       select dimensions in which to select regions of interest
% - 'view&ROI'  select dimensions to view and select ROIs in these
%               dimensions; see 'view' for value format
% - 'display mode'      'time courses' or 'image' (NB: it is possible to
%               specify directly 'time courses' or 'image' without using
%               the 'display mode' argument)
% - 'colormap'  n*3 array or the name of a recognized color map
%               (this automatically set display mode to 'image')
% - 'clipping'  set auto-clipping mode, see 
% - 'controls'  'on'/True (default) or 'off'/False - show/hide the control panel
% ---
% type 'xplor demo' to select a range of demos
% type 'xplor animation' to launch a small animation that illustrates the concepts behind xplor
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
    % xplr.wizard;
    help xplor
    return
elseif nargin == 1 && ischar(varargin{1}) 
    flag = varargin{1};
    switch flag
        case 'demo'
            xplor.demo
            return
        case 'animation'
            xplor.animation
            return
    end
end

% % drag & drop window
% if nargin==0
%     xplor_window()
%     return
% end

% Gather options inside a structure
data = [];
options = struct();
done_with_data = false;
idx = 0;
while idx < length(varargin)
    idx = idx + 1;
    argi = varargin{idx};

    % Data
    if done_with_data
        % argument is not about data but about options
    elseif ischar(argi) || isstring(argi)
        % argument is now a string, i.e. an option name -> wrap-up data
        done_with_data = true;
    elseif isempty(data)
        data = argi;
        continue
    elseif isnumeric(argi)
        % several data arrays give - we will try to concatenate them, for
        % the moment lets just gather them
        if iscell(data)
            data{end+1} = argi;
        elseif isnumeric(data)
            data = {data, argi};
        else
            error argument
        end
        continue
    else
        error argument
    end
    
    % Options
    name = argi;
    if ismember(name, {'time courses', 'image'})
        value = name;
        name = 'display mode';
    else
        if strcmp(name, 'displaymode')
            name = 'display mode';
        end
        idx = idx + 1;
        value = varargin{idx};
    end
    name = brick.strrep(name, ' ', '_', '&', '_and_');
    options.(name) = value;
end

% Concatenate several arrays into a single one
if iscell(data)
    n_data = length(data);
    data_sizes = zeros(n_data, 1);
    for k = 1:n_data
        s = size(data{k});
        if length(s)==2 && s(2)==1, s = s(1); end
        for i = 1:length(s)
            data_sizes(k,i) = s(i);
        end
    end
    same_size = ~any(diff(data_sizes));
    if all(same_size)
        % concatenate in a new dimension
        nd = size(data_sizes,2);
        data = cat(nd+1, data{:});
    elseif sum(~same_size)==1
        % concatenate in the dimension where sizes differ
        d = find(~same_size);
        data = cat(d, data{:});
    else
        % sizes differ in too many dimensions
        error 'Multiple arrays cannot be concatenated because of size mismatch.'
    end
end


% Convert Matlab array to xplr.XData, create headers and launch view
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
    if isvector(data), data = data(:); end % change row vector to column vector
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
V = xplr.View(data, options);

 
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
