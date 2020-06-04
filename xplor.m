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
    %     % flow dataset
    %     [~, ~, ~, dat] = flow;
    %     dat = permute(dat, [1, 3, 2]);
    %     s = size(dat);
    %     head = xplr.Header( ...
    %         {'x', 'px', s(1)}, {'y', 'px', s(2)}, ...
    %         {'axis', num2cell(['a':'y' 'A':'Y'])} ...
    %         );
    %     data = xplr.XData(dat, head, 'Flow Data');
    
    % xplor logo
    logo = fn_readimg(fullfile(fileparts(which('xplor')),'demo','XPLOR logo.png'));
    logo = single(logo)/255;
    logo_rgb = num2cell(logo, [1 2]);
    
    % coordinates in original image: place white sides outside
    white = all(logo == 1, 3);
    [white_x, white_y] = deal(all(white, 2), all(white, 1));
    i_start = find(~white_x, 1, 'first');
    i_end = find(~white_x, 1, 'last');
    j_start = find(~white_y, 1, 'first');
    j_end = find(~white_y, 1, 'last');
    if i_end - i_start > j_end - j_start
        [j_start, j_end] = dealc((j_start+j_end)/2 + [-1 1] * (i_end-i_start)/2);
    else
        [i_start, i_end] = dealc((i_start+i_end)/2 + [-1 1] * (j_end-j_start)/2);
    end
    [xx1, yy1] = ndgrid( ...
        interp1([i_start, i_end], [-1, 1], 1:size(logo,1), 'linear', 'extrap'), ...
        interp1([j_start, j_end], [-1, 1], 1:size(logo,2), 'linear', 'extrap'));
    
    % downsample and rotate
    n2 = 200;
    % coordinates in new image
    [xx2, yy2] = ndgrid(linspace(-1, 1, n2), linspace(-1, 1, n2)); 
    xy2 = [row(xx2); row(yy2)];
    theta = column(linspace(0, 2*pi, n2));
    % convert to coordinates in the original image where to interpolate
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    xy2 = reshape(R * xy2, [n2, 2, n2, n2]);
    xx2 = permute(xy2(:, 1, :, :), [3, 4, 2, 1]);
    yy2 = permute(xy2(:, 2, :, :), [3, 4, 2, 1]);
    % perform interpolation
    dat = zeros(n2, n2, 3, n2);
    for k = 1:3
        dat(:, :, k, :) = interpn(xx1(:, 1), yy1(1, :), logo_rgb{k}, xx2, yy2, 'linear');
    end
    dat(isnan(dat)) = 1;
    
    % XData object
    data = xplr.XData(dat, ...
        {{'x', 'px'}, {'y', 'px'}, {'color', {'r', 'g', 'b'}}, {'z', 'px'}}, ...
        'LOGO');
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
elseif ~isa(data, 'xplr.XData')
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
