function b = reallydlg(varargin)
%REALLYDLG Ask for confirmtaion
%---
% function b = reallydlg(line1,line2,..)
%--
% This is a shortcut for using Matlab questdlg function.
% returns true if 'Yes' has been answered
% argument can be one or several string or cell array of strings
% 
% See also brick.dialog_questandmem

% Thomas Deneux
% Copyright 2009-2017

if nargin==0, help brick.reallydlg, return, end

for i=1:nargin
    if ~iscell(varargin{i}), varargin{i}={varargin{i}}; end
end
question = [varargin{:}];
answer = questdlg(question,'warning','Yes','No','No');
b = strcmp(answer,'Yes');
