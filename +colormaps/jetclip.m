function cm = jetclip(n)
% function cm = jetclip(n)

if nargin<1, n = 10; end
n = min(n,126);
m = floor(126/n);

cm = jet(256); range = [11 246];
cm = interp1(range(1):range(2),cm(range(1):range(2),:),linspace(range(1),range(2),n));
cm = kron(cm,ones(m,1));
cm = [zeros(1,3); cm; ones(1,3)];