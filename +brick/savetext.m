function savetext(a,filename)
% function savetext(a[,filename])
%---
% Input:
% - a           char array
% - filename    file name
%
% See also brick.readtext, brick.readxml, brick.readasciimatrix, brick.readbin

% Thomas Deneux
% Copyright 2005-2017

if nargin==0, help brick.savetext, return, end

if nargin<2, filename=brick.savefile; end

if filename==0, disp(['Could not open file ' filename]), return, end 

fid=fopen(filename,'w');

if ~iscell(a), a = cellstr(a); end
a = brick.strrep(a,'%','%%','\','\\');
a = a(:)';
[a{2,:}] = deal('\n');
s = char([a{:}]);
fprintf(fid,s);

fclose(fid);