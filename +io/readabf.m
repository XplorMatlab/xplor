function x = readabf(file_name)

% file name
if nargin < 1
    file_name = brick.getfile;
end
if ~iscell(file_name)
    file_name = {file_name};
end
nf = length(file_name);

% read file
[data, si, info] = abfload(file_name{1});
if nf > 1
    data(:, :, nf) = 0;
    for i = 2:nf
        data(:, :, i) = abfload(file_name{i});
    end
end

% time info
channel_names = info.recChNames;
if ~strcmp(channel_names, 'Time')
    error 'First channel is not Time, don''t know how to find time information')
end
time = data(:, 1);
data(:, 1) = [];
channel_names(1) = [];
dt = mean(diff(time));

% convert to XPLOR format
if nf ==1 
    x = xplr.XData(data, ...
        {{'time' 's' dt}, ...
        {'channel' channel_names}});
else
    x = xplr.XData(data, ...
        {{'time' 's' dt}, ...
        {'channel' channel_names}, ...
        {'file' brick.fileparts(file_name, 'base')}});
end