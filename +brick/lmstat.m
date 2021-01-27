% Get information about how many Matlab floating licenses are used on the network 

% Thomas Deneux
% Copyright 2013-2021

switch lower(computer)
    case {'pcwin' 'pcwin64'}
        if strcmpi(computer, 'pcwin64')
            w = 'win64';
        else
            w = 'win32';
        end
        attempts = {fullfile(matlabroot, 'bin', w), ...
            fullfile(matlabroot, 'etc', w), ...
            fileparts(matlab.desktop.editor.getActiveFilename)};
        for k=1:length(attempts)
            lmdir = attempts{k};
            ok = exist(fullfile(lmdir,'lmutil.exe'),'file');
            if ok, break, end
        end
        if ~ok, error 'cannot locate file lmutil.exe', end
    otherwise
        lmdir= fullfile(matlabroot, 'etc', lower(computer));
        cmd = ['./' cmd];
end

system(['"' fullfile(lmdir, 'lmutil') '" lmstat -a -c "' fullfile(matlabroot, 'licenses', 'network.lic') '"']);
