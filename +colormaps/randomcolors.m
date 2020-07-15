function col = randomcolors(n,satrange,lumrange)
% function col = randomcolors(n[,saturationrange,luminancerange])
%---
% a colormap with random colors (by default with high saturation and
% luminance)

% input
if nargin<3
    if nargin==2, error 'wrong number of arguments', end
    satrange = [.5 1];
    lumrange = [.9 1];
else
    if isscalar(satrange), satrange = satrange([1 1]); end
    if isscalar(lumrange), lumrange = lumrange([1 1]); end
end

% (random) color parameters
rng('default')
par = rand(3,n)'; % the random numbers are generated in the following order:
% first the 3 parameters of 1st color, then 3 parameters of 2nd color, etc.
% This will ensure that the same initial colors are generated for different
% values of n.

hue = par(:,1);
sat = satrange(1)+diff(satrange)*par(:,2);
lum = lumrange(1)+diff(lumrange)*par(:,3);

% build colors
col = matrix(hsv2rgb(hue,sat,lum));
