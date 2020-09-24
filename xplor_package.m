function xplor_package

%% Create folders

disp 'create build folder'
base_folder = fileparts(which('xplor'));
build_folder = fullfile(base_folder, 'build');
if exist(build_folder,'dir')
    try
        rmdir(build_folder,'s')
    catch ME
        warning(ME.message)
    end
end
mkdir(build_folder)

% disp 'create distrib folder'
% base_folder = fileparts(which('xplor'));
% distrib_folder = fullfile(base_folder, 'distrib');
% if exist(distrib_folder,'dir')
%     try
%         rmdir(distrib_folder,'s')
%     catch ME
%         warning(ME.message)
%     end
% end
% mkdir(distrib_folder)

%% Copy m-files to the build folder while flatening the sub-directory and namespace structure
send_m_file(fullfile(base_folder, 'xplor.m'), build_folder)
copyfile(fullfile(base_folder, 'xplor parameters.xml'), build_folder)
send_folder(fullfile(base_folder, '+brick'), build_folder, 'brick_')
copyfile(fullfile(base_folder, '+brick', '*.png'), build_folder)
send_folder(fullfile(base_folder, '+brick', 'private'), build_folder, '')
send_folder(fullfile(base_folder, '+colormaps'), build_folder, 'colormaps_')
send_folder(fullfile(base_folder, '+io'), build_folder, 'io_')
send_folder(fullfile(base_folder, '+xplor'), build_folder, '')
send_folder(fullfile(base_folder, '+xplor', '+demo'), build_folder, 'demo_')
send_folder(fullfile(base_folder, '+xplr'), build_folder, '')
send_folder(fullfile(base_folder, '+xplr', 'private'), build_folder, '')

%% Compile
file_list = dir(build_folder);


%% Functions
function send_m_file(file_name, build_folder, prefix)

% display
f = brick.fileparts(file_name, 'name');
disp(f)

% read
str = fileread(file_name);

% replace
str = regexprep(str, '(xplr|xplor)\.', '');
str = regexprep(str, '(demo|brick|colormaps|io)\.', '$1_');

% write
if nargin >= 3
    fun = brick.fileparts(file_name, 'base');
    if strcmp(str(1:8), 'classdef')
        str = regexprep(str, ['(?<![a-zA-Z_0-9])' fun '(?![a-zA-Z_0-9])'], [prefix fun]);
    elseif strcmp(str(1:8), 'function')
    end
    new_file_name = fullfile(build_folder, [prefix f]);
else
    new_file_name = fullfile(build_folder, f);
end
fid = fopen(new_file_name, 'w');
fwrite(fid, str);
fclose(fid);

function send_folder(folder, build_folder, prefix)

d = dir(fullfile(folder, '*.m'));
for i = 1:length(d)
    send_m_file(fullfile(folder, d(i).name), build_folder, prefix)
end






