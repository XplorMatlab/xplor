function varargout = regexptokens(a,expr,ktoken)
%REGEXPTOKENS Get the tokens of a regexp as a simple cell array
%---
% function tokens = regexptokens(a,expr)
% function [tok1 ... tokn] = regexptokens(a,expr)
% function tokk = regexptokens(a,expr,ktoken)
%---
% a wrapper of regexp that returns the tokens in a string where expression
% 'expr' is assumed to occur exactly one time
% returns the unique token if there is only one, or a cell array of tokens
% if there are several

% Thomas Deneux
% Copyright 2015-2017

if ~iscell(a)
    a = {a};
end
n = numel(a);
if n==0
    varargout = cell(1, max(nargout,1));
    [varargout{:}] = deal(cell(1,0));
    return
end
tokens = regexp(a,expr,'tokens');       % cell array (n) of cell arrays (#match) of cell arrays (#tokens) of strings
if isempty(tokens) || isempty(tokens{1})
    error 'no match'
end
ntoken = length(tokens{1}{1});
if ntoken==0                
    disp 'brick.regexptokens: no token, returning full match'
    tokens = brick.map(@(x){x}, regexp(a,expr,'match'), 'cell');
    ntoken = 1;
end
% out is a cell array (#tokens) of cell arrays (n) of strings
out = cell(1, ntoken);
[out{:}] = deal(cell(1,n));
for i=1:n
    tokens_i = tokens{i};               % one or zero element depending on match or not
    if length(tokens_i)>1
        error 'more than one match'
    elseif isempty(tokens_i)
        error 'no match'
    elseif length(tokens_i{1}) ~= ntoken
        error 'not the same number of tokens in every string, can this really happen?!'
    end
    tokens_i = tokens_i{1};
    for k=1:ntoken
        out{k}{i} = tokens_i{k}; 
    end
end

% Simplify output if n==1 or ntoken==1 or specific token requested
if n==1
    out = brick.map(@(x)x{1}, out, 'cell');    % cell array (#tokens) of strings
end
if isempty(out)
    % there was not a single match, so we don't know what was the supposed
    % number of tokens, get it from the number of requested outputs or
    % assume 1
    ntoken = max(nargout,1);
    out = cell(1, ntoken);
    [out{:}] = deal(cell(1,n));
end
if ntoken==1
    out = out{1};                   % string (if n==1) or cell array (n) of strings
elseif nargin>=3
    out = out{ktoken};
end

% Output
if nargout <= 1
    varargout = {out};
else
    if nargout~=length(out), error 'number of outputs does not match number of tokens', end
    varargout = out;
end
