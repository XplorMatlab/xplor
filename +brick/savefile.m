function [filename filterindex] = savefile(varargin)
%SAVEFILE User select file for saving and remember last containing folder
%---
% function [filename filterindex] = savefile([filter[,title]])
%--
% synonyme de "[filename filterindex] = brick.getfile('SAVE',[filter[,title]])"
% 
% See also brick.getfile

% Thomas Deneux
% Copyright 2003-2017

[filename filterindex] = brick.getfile('SAVE',varargin{:});
