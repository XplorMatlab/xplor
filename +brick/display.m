function display(varargin)
% function display(x[,format])
% function display(varname,x[,format])
%---
% Display a variable name and value in a single line

% Thomas Deneux
% Copyright 2015-2017

if nargin==0, help brick.display, return, end

% Input
xname = []; x = []; format = [];
for i=1:nargin
    a = varargin{i};
    if isnumeric(a)
        x = a;
    elseif any(a=='%')
        format = a;
    elseif isempty(xname)
        xname = a;
    elseif isempty(x)
        x = evalin('caller',a);
    end
end
if isempty(xname)
    xname = inputname(1);
elseif isempty(x)
    x = evalin('caller',xname);
end

% Display
if isempty(format)
    disp([xname ': ' num2str(x)])
else
    disp([xname ': ' num2str(x,format)])
end

       
