function log_to_file(str, command)
% function log_to_file(str)
% function log_to_file([], 'init')
%---
% log to file 'xplor.log' in the main xplor folder

file_name = fullfile(fileparts(which('xplor')), 'xplor.log');

if isempty(str)
    switch command
        case 'init'
            if exist(file_name, 'file')
                delete(file_name)
            end
            str = 'init log';
        otherwise
            error('unknown command ''%s''', command)
    end
end

% write to log file
fid = fopen(file_name, 'a');
fprintf(fid, [str '\n']);
fclose(fid);
