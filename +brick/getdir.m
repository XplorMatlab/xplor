function filename = getdir(varargin)
%GETDIR Select directory and remember last containing folder
%---
% function dirname = getdir(title[,dirname])
%--
% synonyme de "filename = brick.getfile('DIR',title)"

% Thomas Deneux
% Copyright 2003-2017

filename = brick.getfile('DIR',varargin{:});
