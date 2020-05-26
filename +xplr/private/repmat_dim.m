function y = repmat_dim(x,dim,rep)
% function y = repmat_dim(x,dim,rep)
%---
% Replicate and tile arrays as Matlab's repmat, but only in the specified
% dimension(s)
% 
% See also repmat

if isempty(dim)
    y = x;
    return
end
reps = ones(1,max(max(dim),2));
reps(dim) = rep;
y = repmat(x, reps);