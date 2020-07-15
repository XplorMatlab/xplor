function y = nrms(x,dim)
% function y = nrms(x,dim)
%---
% returns RMS while ignoring NaNs
%
% see also nmean

% Thomas Deneux
% Copyright 2015-2017


if nargin==1
  y = sqrt(brick.nmean(x .* conj(x)));
else
  y = sqrt(brick.nmean(x .* conj(x), dim));
end

