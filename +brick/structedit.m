function s_out = structedit(s,varargin)
%STRUCTEDIT Edit a structure with a GUI brick.interface; this is a wrapper of brick.control
%---
% function s = structedit(s[,spec][,hp][,other control arguments...][,'memorize'])
% function s = structedit(fieldname1,value1,fieldname2,value2,...)
%---
% Allow user to interactively modify a structure, and returns the modified
% structure; returns an empty array if the figure is closed
% this function is a wrapping of control.
% 'memorize' option will cause the result to be saved on the disk so that
% if the function is called again with the same default value, this default
% value will be replaced by the lastly saved value.
% 
% See all brick.control, brick.input, brick.reallydlg

% Thomas Deneux
% Copyright 2007-2017

if ischar(s)
    % build a 3-elements structure
    structdef = [s varargin];
    for k=2:2:length(structdef)
        a = structdef{k};
        if iscell(a)
            [a{end+1:3}] = deal([]);
        else
            a = {a [] []}; %#ok<AGROW>
        end
        structdef{k} = a;
    end
    s = struct(structdef{:});
    varargin = {};
end

% replace default with previously memorized value?
do_memorize = ~isempty(varargin) && strcmp(varargin{end}, 'memorize');
if do_memorize
    varargin(end) = []; 
    mem_hash = ['structedit_' brick.hash(s(1))];
    s_saved = brick.userconfig(mem_hash);
    if ~isempty(s_saved)
        s(1) = s_saved;
    end
end
X = brick.control(s,varargin{:},'ok');
s_out = [];

% nested function: update of s
function update
    s_out = X.s;
end

addlistener(X,'OK',@(u,e)update());
update()

% wait that the figure will be destroyed
waitfor(X.hp)
drawnow % otherwise figure will stay on if heavy computation continue afterwards

% memorize value?
if do_memorize
    brick.userconfig(mem_hash, s_out)    
end

end
        
        

