function b = same_size(sz1,sz2)
% function b = same_size(sz1,sz2)
%---
% compare two size vectors (be careful, the inputs should not be arrays
% whose size we want to chek, but the size vectors; otherwise use check_size)
%
% See also check_size

n = max(length(sz1),length(sz2));
sz1(end+1:n) = 1;
sz2(end+1:n) = 1;
b = all(brick.row(sz1)==brick.row(sz2));
