function xplor_package

%% Create folders

disp 'create build folder'
base_folder = fileparts(which('xplor_package'));
build_folder = fullfile(base_folder, 'build');
if exist(build_folder,'dir')
    try
        rmdir(build_folder,'s')
    catch ME
        warning(ME.message)
    end
end
mkdir(build_folder)


%% Copy m-files to the build folder while flatening the sub-directory and namespace structure
disp 'copy files'
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
copyfile(fullfile(base_folder, '+xplr', 'wizard.html'), build_folder)
send_folder(fullfile(base_folder, '+xplr', 'private'), build_folder, '')
send_m_file(fullfile(base_folder, '..', 'MatlabPerso', 'work', 'optimage_checklicense.m'), build_folder)

%% Compile (manually or using the Application Compiler app)

DOMANUAL = false;

if DOMANUAL
    
    %% Get dependencies
    disp 'get dependencies'

    cd(build_folder)
    m_files = matlab.codetools.requiredFilesAndProducts('xplor.m');
    idx = brick.find(strfind(m_files, 'xplor.m'));
    m_files = m_files([idx, setdiff(1:length(m_files), idx)]);

    file_list = dir(build_folder);
    file_list([file_list.isdir]) = [];
    file_list = {file_list.name};
    is_m_file = ~brick.isemptyc(regexp(file_list, '.m$'));
    other_files = file_list(~is_m_file);
    other_files = [repmat({'-a'},[1 length(other_files)]); other_files];

    %% Compile

    disp 'compile'
    cd(build_folder)
    mkdir distrib
    mcc('-e', '-R', '-logfile,xplor.log', m_files{:}, other_files{:})

    %% Make installer

    disp 'create distrib folder'
    distrib_folder = fullfile(base_folder, 'distrib');
    if exist(distrib_folder,'dir')
        try
            rmdir(distrib_folder,'s')
        catch ME
            warning(ME.message)
        end
    end
    mkdir(distrib_folder)
    cd(distrib_folder)

    disp 'make installer'
    install_files = { ...
        fullfile(build_folder, 'xplor.exe'), ...
        fullfile(base_folder, 'demo data')};
    MCR_file = fullfile(build_folder, 'requiredMCRProducts.txt');
    logo_file = fullfile(base_folder, 'XPLOR logo.png');

    cd(distrib_folder)
    compiler.package.installer( ...
        install_files, MCR_file, ...
        'ApplicationName', 'XPLOR',...
        'AuthorName', 'Thomas Deneux',...
        'InstallerName', 'XPLOR Installer',...
        'Summary', 'XPLOR your data!', ...
        'InstallerIcon', logo_file, ...
        'InstallerLogo', logo_file, ...
        'InstallerSplash', logo_file, ...
        'RunTimeDelivery', 'web', ...
        'Shortcut', fullfile(build_folder, 'xplor.exe'))
    
else
    %% use the Application Compiler app
    disp 'Application Compiler app'
    cd(build_folder)
    applicationCompiler(fullfile(base_folder, 'xplor.prj'))

end

%% Functions
function send_m_file(file_name, build_folder, prefix)

% display
source_folder = brick.fileparts(fileparts(file_name), 'name');
f = brick.fileparts(file_name, 'name');
% disp(f)

% read
str = fileread(file_name);

% replace
str = regexprep(str, '(xplr|xplor)\.', '');
str = regexprep(str, '(demo|brick|colormaps|io)\.', '$1_');
if ismember(source_folder, {'+demo' '+brick' '+colormaps' '+io'})
    % replace function name in the first line
    base = brick.fileparts(f, 'base');
    prefix = [source_folder(2:end) '_'];
    str = regexprep(str, ['(?<=^[^\n]*)' base], [prefix base]);
elseif strcmp(f, 'wizard.m')
    str = regexprep(str, '\+xplr/', '');
end

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





