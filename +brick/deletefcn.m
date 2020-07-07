function deletefcn(hu,deletefcn)
% function deletefcn(hu,deletefcn)
%---
% Set the 'DeleteFcn' property of a graphic object in such a way that
% several functions can be executed upon its deletion.

% Thomas Deneux
% Copyright 2015-2017

warning 'function deletecfn(hu,deletefcn) is deprecated, use addlistener(hu,'ObjectBeingDestroyed',deletefcn); instead'

% Multiple objects
if iscell(hu) || ~isscalar(hu)
    if ~iscell(hu), hu = num2cell(hu); end
    for i=1:numel(hu)
        brick.deletefcn(hu{i},deletefcn)
    end
    return
end

% Listen to object deletion
addlistener(hu,'ObjectBeingDestroyed',deletefcn);
    
