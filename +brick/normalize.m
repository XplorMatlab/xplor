function x = normalize(x,dim,flag)
%NORMALIZE Normalize data based on averages in some specific dimension
%---
% function x = normalize(x,dim,flag)
%---
% Input:
% - x       array
% - dim     dimensions on which to operate; can be a cell array for
%           multiple actions
% - action  'div'[default], 'sub', 'std' or 'zscore', 'norm2', 'detrend', 
%           'proba' (divide by the sum rather than by the mean), 'max',
%           'minmax'
%           can be a cell array for multiple actions 

% Thomas Deneux
% Copyright 2005-2017

if nargin<2, dim=1; end
if nargin<3, flag='div'; end

if isempty(x), return, end

if iscell(dim)
    nrep = length(dim);
    if ~iscell(flag), [tmp{1:nrep}] = deal(flag); flag=tmp; end
    for i=1:nrep
        x = brick.normalize(x,dim{i},flag{i});
    end
    return
end
    
nd = length(dim);
if ~isa(x,'single'), x = double(x); end
m = x; M = x;
for k=1:nd
    switch flag
        case 'max'
            m = max(m,[],dim(k));
        case 'proba'
            m = sum(m,dim(k)); 
        case 'minmax'
            m = min(m,[],dim(k));
            M = max(M,[],dim(k));
        otherwise
            m = brick.nmean(m,dim(k));
    end
end
switch flag
    case {'div','/','proba','max'}
        x = brick.div(x,m);
    case 'minmax'
        x = brick.div(brick.subtract(x,m),M-m);
    case {'sub','-'}
        x = brick.subtract(x,m);
    case {'std','zscore'}
        x = brick.subtract(x,m);
        m2 = x.^2;
        for k=1:nd
            m2 = brick.nmean(m2,dim(k));
        end
        x = brick.div(x,sqrt(m2));
    case 'norm2'
        m2 = x.^2;
        for k=1:nd
            m2 = sum(m2,dim(k));
        end
        x = brick.div(x,sqrt(m2));
    case 'detrend'
        if nd~=1, error('detrend should be along only one dimension'), end
        nt = size(x,dim);
        tt = shiftdim((1:nt)',1-dim); tt = tt-mean(tt); tt = tt/sqrt(sum(tt.^2));
        a = sum(brick.mult(x,tt),dim);
        x = x - brick.add(brick.mult(a,tt),m);
    otherwise
        error(['unknown flag ''' flag ''''])
end

