function savexml(fname,s)
% function savexml(fname,s)
%---
% This function is a wrapper for function xml_write.m downloaded on Matlab
% File Exchange

% Jarek Tuszynski (function xml_write.m)
% Copyright 2007-2012
% Thomas Deneux
% Copyright 2007-2017

ext = brick.fileparts(fname,'ext');
if isempty(ext), fname = [fname '.xml']; end
xml_write(fname,s);
