function s = num2str(x,varargin)
%NUM2STR Convert numeric to char, unless input is already char!, can return a cell array
%---
% function s = num2str(x[,'cell'][,'quotestrings'][,'format'])
%---
% convert numerical value x into a string representation... unless s is
% already a character array!
% 
% use the 'cell' flag to make the output a cell array (the same size as x)
%
% See also brick.str2double, brick.strcat, brick.chardisplay, brick.idx2str

% Thomas Deneux
% Copyright 2007-2017

if nargin==0, help brick.num2str, return, end

% flag
docell = iscell(x); format = {}; doquotestrings = false;
for k=1:length(varargin)
    a = varargin{k};
    if strcmp(a,'cell')
        docell = true;
    elseif strcmp(a,'quotestrings')
        doquotestrings = true;
    else
        format = {a};
    end
end

% make a cell array first
if iscell(x)
    % nothing to do
elseif ischar(x)
    x = {x};
elseif isnumeric(x) || islogical(x) || isdatetime(x) || isduration(x)
    x = num2cell(x);
else
    error argument
end

% conversion
s = x;
for i=1:numel(x)
    xi = x{i};
    if isnumeric(xi) || islogical(xi)
        s{i} = num2str(xi, format{:});
    elseif isdatetime(xi) || isduration(xi)
        if ~isempty(format)
            xi.Format = format{1};
        end
        s{i} = char(xi);
    elseif doquotestrings
        s{i} = ['"' xi '"'];
    end
end

% convert back to char if required
if ~docell
    % add spaces in the 2nd dimension
    if size(s,2)>1
        s = brick.interleave(2,s,repmat({' '},size(s)));
        s(:,end,:) = [];
    end
    s = cell2mat(s);
end
