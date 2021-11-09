function [filename filterindex] = getfile(varargin)
%GETFILE Select file and remember the containing folder of the last selected file
%---
% function [filename filterindex] = getfile(['READ|SAVE',][filter[,title]])
% function filename = getfile('DIR',title)
% function rep = getfile('REP')
% function getfile('REP',rep)
%---
% returns file name with full path
% brick.getfile('REP') returns the current directory
% brick.getfile('REP',rep) sets the current directory
% as a utility, brick.getfile('cd') is going to the current directory
% 
% See also brick.savefile

% Thomas Deneux
% Copyright 2003-2017

persistent rep

% cd to persistent directory
swd=pwd;
if isempty(rep), rep = pwd; end
try cd(rep), catch, rep = pwd; end

% Input
arg = varargin;
if (length(arg)>=1 && ischar(arg{1}) && brick.ismemberstr(arg{1},{'READ','SAVE','DIR','REP','cd'}))
    flag=arg{1}; 
    arg={arg{2:end}}; 
else
    flag='READ';
end
arg0 = arg;
if length(arg)<1 || isempty(arg{1}), arg{1}='*.*'; end
if length(arg)<2 || isempty(arg{2}), arg{2}='Select file'; end

% cd to provided path?
if brick.ismemberstr(flag,{'READ' 'SAVE'})
    fil = arg{1};
    if iscell(fil), fil = fil{1}; end
    p = fileparts(fil);
    if ~isempty(p)
        cd(p)
        fil = brick.fileparts(fil,'name');
        if iscell(arg{1}), arg{1}{1}=file; else arg{1}=fil; end
    end
end

% Matlab running in no display mode?
if ~feature('ShowFigureWindows') && ismember(flag,{'READ' 'SAVE' 'DIR'})
    disp 'brick.getfile: cannot prompt user for file or directory because Matlab is in no display mode'
    filename = 0;
    filterindex = 0;
    return
end

% different actions
switch flag
    case 'READ'
        if ischar(arg{1}) && ~any(arg{1}=='*'), arg{1} = {arg{1}, ['*' brick.fileparts(arg{1},'ext')]}; end
        if any(strfind(computer, 'MAC')) && length(arg)>=2
            % Displqy prompt because it will not appear in the getfile
            % window for MAC
            disp(arg{2})
        end
        [filename pathname filterindex] = uigetfile(arg{:},'MultiSelect','on');
        if iscell(filename), filename = char(filename{:}); end 
        if filename
            filename = [repmat(pathname,size(filename,1),1) filename];
            rep = pathname;
        end
    case 'SAVE'
        if any(strfind(computer, 'MAC')) && length(arg)>=2
            % Displqy prompt because it will not appear in the getfile
            % window for MAC
            disp(arg{2})
        end
        [filename pathname filterindex] = uiputfile(arg{:});
        if filename
            filename = fullfile(pathname,filename);
            rep = pathname;
        end
    case 'DIR'
        switch length(arg0)
            case 0
                arg = {[] 'Select directory'};
            case 1
                arg = [{[]} arg0];
            case 2
                arg = arg0([2 1]);
            otherwise
                error('too many arguments')
        end
        if any(strfind(computer, 'MAC')) && length(arg)>=2
            % Displqy prompt because it will not appear in the getfile
            % window for MAC
            disp(arg{2})
        end
        filename = uigetdir(arg{:});
        if filename
            rep = fileparts(filename);
        end
    case 'REP'        
        % special case: access/set the persistent directory name rep
        if nargin==1
            filename = rep;
        elseif nargin==2
            rep = varargin{2};
        end
    case 'cd'
end

% cd back to initial directory
if ~strcmpi(flag,'cd'), cd(swd), end
