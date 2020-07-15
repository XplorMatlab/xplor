function hl = drawpoly(poly,varargin)
%DRAWPOLY Shortcut for line(poly(:,1),poly(:,2))
%---
% function hl = drawpoly(poly[,'close'],line options...)
%---
% shortcut for line(poly(1,:),poly(2,:),varargin{:})
%
% See also brick.circle

% Thomas Deneux
% Copyright 2006-2017

if nargin==0, help brick.drawpoly, end

if isa(poly,'selectionND')
    poly = poly.convert('poly2D');
    poly = poly.poly.points;
end

if size(poly,1)~=2
    if size(poly,2)==2
        poly = poly';
    else
        error('poly should have two rows')
    end
end
doclose = ~isempty(varargin) && strcmp(varargin{1},'close');
if doclose, varargin(1) = []; end
opt = brick.linespecs(varargin{:});
if doclose
    hl = line(poly(1,[1:end 1]),poly(2,[1:end 1]),opt{:});
else
    hl = line(poly(1,:),poly(2,:),opt{:});
end
if nargout==0, clear hl, end
