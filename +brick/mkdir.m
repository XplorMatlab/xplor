function mkdir(D)
%MKDIR Create a directory if it does not exist
%---
% function mkdir(D)
%---
% creates directory D if it does not exist

% Thomas Deneux
% Copyright 2004-2017

if nargin==0, help brick.mkdir, end

if ~exist(D,'dir'), mkdir(D), end
