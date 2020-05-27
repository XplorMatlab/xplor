function b = check_size(x,sz)
% function b = check_size(x,sz)
%---
% check that array has the desired size
%
% See also check_size

b = same_size(size(x),sz);

if nargout==0 
    if b
        clear b
    else
        error(['Array has size ' num2str(strictsize(x),'%i ') ', expected ' num2str(sz,'%i ')])
    end
end