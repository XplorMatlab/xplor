function fname = autofigname(hf)
% function fname = autofigname(hf)
% function autofigname('setfolder')
% ---
% creates an automatic name for file where to save figure hf (without
% extension however)

% Thomas Deneux
% Copyright 2012-2017

if nargin==0, help brick.autofigname, return, end

persistent savedir
if isempty(savedir) || strcmp(hf,'setfolder')
    d = brick.userconfig('brick.autofigname');
    d = uigetdir(d,'Select folder where to save figures');
    if ~ischar(d)
        error('No folder was selected!')
    elseif ~exist(d,'dir')
        error('''%s'' is not a valid folder name!',d)
    end
    brick.userconfig('brick.autofigname',d)
    savedir = d;
    if strcmp(hf,'setfolder'), return, end
end

figname = get(hf,'name');
hfnum=get(hf,'Number');
if isempty(figname), figname = ['Figure' num2str(hfnum)]; end
fname = fullfile(savedir, ... 
    [regexprep(figname,'[ |:/\\\*\+<>]','_')  datestr(now,'_yymmdd_HHMMSS')]);
