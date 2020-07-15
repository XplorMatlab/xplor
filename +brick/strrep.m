function str = strrep(str,varargin)
%STRREP Replace several text sequences in string
%---
% function str = strrep(str,pattern1,rep1,pattern2,rep2,...)
%---
% Same as Matlab strrep, but allows multiple pattern replacements

% Thomas Deneux
% Copyright 2015-2017

if nargin==0, help brick.strrep, return, end

if ~mod(nargin,2), error 'pattern and replacement strings must come as pairs', end

for k=1:2:nargin-2, str = strrep(str,varargin{[k k+1]}); end
