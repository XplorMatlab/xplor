function y = subsref_dim(x, dim, ind)
% function y = subsref_dim(x, dim, ind)
%---
% Simplified version of Matlab's subsref to access sub-portion of array in
% the specified dimension(s);
%
% Example:
% x = [1 2 3; 4 5 6]
% 
% x =
% 
%      1     2     3
%      4     5     6
% 
% y = subsref_dim(x,2,[1 3])
% 
% y =
% 
%      1     3
%      4     6
%
% See also subsref, subsasgn_dim, repmat_dim

if isempty(dim)
    nd = ndims(x);
else
    nd = max(ndims(x), max(dim));
end
if ~iscell(ind)
    if isscalar(dim)
        ind = {ind};
    elseif length(dim) == length(ind)
        ind = num2cell(ind);
    else
        error 'ind must be an array or cell array of same length as dim'
    end
end
inds = repmat({':'},1,nd);
inds(dim) = ind;
subs = substruct('()',inds);
y = subsref(x, subs);
