function cm = bluered(n)

if nargin < 1, n = size(get(gcf,'colormap'),1); end

k = ceil(n/2);

blue = [linspace(0,1,k)'*[1 1] ones(k,1)];
red = [ones(k,1) linspace(1,0,k)'*[1 1]];

if n==2*k
    cm = [blue; red];
elseif n==2*k-1
    cm = [blue; red(2:end,:)];
end
    
cm = exp(10*(cm-.5))./(1+exp(10*(cm-.5)));