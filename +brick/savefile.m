function [filename, filterindex] = savefile(varargin)
%SAVEFILE User select file for saving and remember last containing folder
%---
% function [filename filterindex] = savefile([filter[,title]])
%--
% synonyme de "[filename filterindex] = brick.getfile('SAVE',[filter[,title]])"
% 
% See also brick.getfile

% Thomas Deneux
% Copyright 2003-2017

[filename, filterindex] = brick.getfile('SAVE',varargin{:});

% Update file name according to filter, i.e. add initial and final part
% (extension) if needed.
if ~isequal(filename, 0) && length(varargin) >= 1
    filter = varargin{1};
    filter = brick.fileparts(filter, 'name'); % remove folder component if any
    if iscell(filter), error 'not implemented yet', end
    filter = char(filter);
    idx_star = find(filter == '*');
    if ~isempty(idx_star)
        bef = filter(1:idx_star-1);
        aft = filter(idx_star+1:end);
        [folder, base] = brick.fileparts(filename, 'path', 'name');
        if ~isempty(bef) && ~any(regexp(base, ['^' bef]))
            base = [bef base];
        end
        if ~isempty(aft) && ~any(regexp(base, [aft '$']))
            base = [base aft];
        end
        filename = fullfile(folder, base);
    end
end
