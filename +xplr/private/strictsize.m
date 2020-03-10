function s = strictsize(x,ndim)

s = size(x);
nd = find(s~=1,1,'last');
if isempty(nd), nd = ~isempty(x); end
s(nd+1:end) = [];
if nargin>=2, s(end+1:ndim)=1; end

end
