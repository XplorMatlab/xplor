function fname = fileext(fname,ext)
% function fname = fileext(fname,ext)
%---
% Set (or replace) extension of file name

% Thomas Deneux
% Copyright 2003-2017

% input
if ext(1)~='.', ext = ['.' ext]; end

% output
fname = [brick.fileparts(fname,'noext') ext];
