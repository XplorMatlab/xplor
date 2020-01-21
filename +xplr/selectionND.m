classdef selectionND < xplr.graphnode
    % function sel = selectionND('type',data[,sizes])
    %---
    % Type can be 'point1D', 'line1D', 'point2D', 'line2D', 'poly2D',
    % 'rect2D', 'ellipse2D', 'ring2D'. Use selectionND('empty',nd) to get
    % an empty selection.
    % If 'sizes' is specified, indices are computed according to this data
    % sizes.
    % 
    % 'selectionND' is a handle class to make a = b faster (no copy), but
    % be careful when using it
    %
    % Syntax for defining selections:
    % - selectionND('point1D',[x1 x2 ...])
    % - selectionND('line1D',[a1 b1 a2 b2 ...])
    % - selectionND('point2D',[x1 x2 ...; y1 y2 ...])
    % - selectionND('line2D',???)
    % - selectionND('poly2D',[x1 x2 ...; y1 y2 ...])
    % - selectionND('rect2D',[x y w h])
    % - selectionND('ellipse2D',{[xc yc],r})            (circle) coordinates of center + radius
    % - selectionND('ellipse2D',{[xc yc],[xu yu],e})    (ellipse) center + principal verctor + eccentricity
    % - selectionND('ring2D',{[xc yc],[xu yu],e,r})     + relative radius
    % - selectionND('ring2D',{[xc yc],[xu yu],[e r]})   
    % - selectionND('all1D|all2D')
    %
    % - selectionND('seltype',{def1 def2 ...})   defines a multiple selection at once
    properties
        active = true;
    end    
    properties (SetAccess='private')
        id      % 2-element vector represent the 'identity' of the selection
                % * all selections in different displays and different
                % referentials which represent the same region share the
                % same id numbers
                % * when the region is changed, the first id number
                % remains the same, but the second is changed
        nd
        poly = struct( ...
            'type',     '', ...
            'points',   {}, ...
            'vectors',  {}, ...
            'logic',    {}, ...
            'special',  {} ...
            );  % empty structure at init
        datasizes = [];
        dataind = 0;
    end
    properties (Dependent, SetAccess='private')
        type % either a possible value for poly.type (e.g. point1D, ellipse2D), or 'mixed' if several types are present
        mask
    end
    
    % Constructor + Load + Display
    methods
        function sel = selectionND(type,data,sizes)
            % function sel = selectionND('type',data[,sizes])
            % function sel = selectionND('empty',nd)
            %---
            % selection ID is made of 2 numbers: the first remains the same as long as
            % the selection is living; the second is changed each time the seleciton is
            % modified
            
            % empty selection (invalid)
            if nargin==0
                sel.nd = 0;
                return
            end
            
            % multiple selection
            if (fn_ismemberstr(type,{'point1D' 'line1D' 'point2D' 'poly2D' 'line2D' 'rect2D'}) && iscell(data)) ...
                    || (fn_ismemberstr(type,{'ellipse2D' 'ring2D'}) && iscell(data{1}))
                if nargin>=3, siz = {sizes}; else siz = {}; end
                nsel = numel(data);
                sel(nsel) = sel; % initialize vect
                for i=1:nsel, sel(i) = selectionND(type,data{i},siz{:}); end
                return
            end
            
            % random ID
            sel.id = rand(1,2);
            
            % number of dimension
            switch type
                case {'point1D' 'line1D' 'all1D'}
                    sel.nd = 1;
                case {'point2D' 'poly2D' 'line2D' 'rect2D' 'ellipse2D' 'ring2D' 'all2D'}
                    sel.nd = 2;
                case 'empty'
                    sel.nd = data;
                    return
                otherwise
                    error('unknown selection type ''%s''',type)
            end
            
            % Syntax for defining selections:
            % - selectionND('point1D',[x1 x2 ...])
            % - selectionND('line1D',[a1 b1 a2 b2 ...])
            % - selectionND('point2D',[x1 x2 ...; y1 y2 ...])
            % - selectionND('line2D',[xa1 xb1 xa2 xb2 ...; ya1 yb1 ya2 yb2 ...])
            % - selectionND('poly2D',[x1 x2 ...; y1 y2 ...])
            % - selectionND('rect2D',[x y w h])
            % - selectionND('ellipse2D',{[xc yc],[xu yu]})      coordinates of center + principal vector
            % - selectionND('ellipse2D',{[xc yc],[xu yu],e})    + eccentricity
            % - selectionND('ring2D',{[xc yc],[xu yu],e,r})     + relative radius
            % - selectionND('ring2D',{[xc yc],[xu yu],[e r]})
            % - selectionND('all1D|all2D')
            
            % set sel
            sel.poly(1).type = type;
            sel.active = true;
            switch type
                case 'point1D'
                    sel.poly.points = data(:)';
                    sel.poly.vectors = zeros(1,0);
                case 'line1D'
                    lines = data(:)';
                    if mod(length(lines),2) || any(diff(lines(:))<=0)
                        error('set of lines should come as ordered non-intersecting segments')
                    end
                    nline = length(lines)/2;
                    lines = reshape(lines,[1 2 nline]);
                    for i=1:nline
                        sel.poly(i).points = lines(:,:,i);
                        sel.poly(i).vectors = zeros(1,0);
                    end
                case 'point2D'
                    if size(data,1)==1; data = data'; end
                    if size(data,1)~=2, error('data should have 2 rows'), end
                    sel.poly.points = data;
                    sel.poly.vectors = zeros(2,0);
                case 'line2D'
                    if size(data,1)==1; data = data'; end
                    if size(data,1)~=2, error('data should have 2 rows'), end
                    if mod(size(data,2),2), error('set of lines improperly defined'), end
                    nline = size(data,2)/2;
                    lines = reshape(data,[2 2 nline]);
                    for i=1:nline
                        sel.poly(i).points = lines(:,:,i);
                        sel.poly(i).vectors = zeros(2,0);
                    end
                case 'poly2D' 
                    if size(data,1)~=2, error('data should have 2 rows'), end
                    if ~all(data(:,1)==data(:,end)), data = data(:,[1:end 1]); end
                    [px py] = deal(data(1,:),data(2,:)); % poly2cw(data(1,:),data(2,:));
                    sel.poly.points = [px(:)'; py(:)'];
                    sel.poly.vectors = zeros(2,0);
                case 'rect2D'
                    if ~isvector(data) || length(data)~=4, error('data should be a 4-element vector (x,y,w,h)'), end
                    if any(data(3:4)<0), error('width and height must be >= 0'), end
                    sel.poly.points = [data(1); data(2)];
                    sel.poly.vectors = [data(3); data(4)];
                case {'ellipse2D' 'ring2D'}
                    if ~iscell(data) 
                        error 'ellipse or ring description must be a cell array';
                    end
                    % center and radius
                    if length(data)<2 || numel(data{1})~=2 || numel(data{2})>2
                        error 'center or radius is ill-defined'
                    end
                    [c, u] = deal(data{1:2});
                    c = c(:); u = u(:);
                    if isscalar(u), u = [u; 0]; end
                    % eccentricity and ring secondary radius
                    switch length(data)
                        case 2
                            if strcmp(type,'ring2D'), error 'eccentricity and secondary radius missing for ring description', end
                            e = 1;
                            logic = e;
                        case 3
                            logic = data{3};
                        case 4
                            [e r] = deal(data{3:4});
                            logic = [e r];
                        otherwise
                            error 'too many elements in shape description'
                    end
                    if length(logic)~=fn_switch(type,'ellipse2D',1,'ring2D',2)
                        error 'wrong shape definition'
                    end
                    % that's it
                    sel.poly.points = c;
                    sel.poly.vectors = u;
                    sel.poly.logic = logic;
                    sel.poly.special = EllipseVector2Sym(u,logic(1));
                case {'all1D' 'all2D'}
                    % nothing need to be set!
                otherwise
                    error('unknown type ''%s''',type)
            end
            
            % indices
            if nargin==3
                sel = ComputeInd(sel,sizes);
            end
        end
        function disp(sel)
            warning('off','MATLAB:structOnObject')
            if isempty(sel)
                disp('empty selectionND object')
            elseif ~isscalar(sel)
                s = size(sel);
                str = cellstr(num2str(s'))';
                [str{2,:}] = deal('x');
                str = [str{1:end-1} ' selectionND object'];
                fprintf('%s\n\nProperties:\n',str)
                F = fieldnames(sel);
                nF = length(F);
                for k=1:nF
                    fprintf('    %s\n',F{k});
                end
            else
                strucdisp(struct(sel))
            end
            warning('on','MATLAB:structOnObject')
        end
    end
    methods (Static)
        function sel = loadobj(sel)
            % previous version might not have a poly.special field
            if ~isfield(sel.poly,'special')
                if isempty(sel.poly)
                    error('programming: a selection cannot have an empty poly')
                end
                sel.poly(1).special = [];
                for k=1:length(sel.poly)
                    a = sel.poly(k);
                    if strcmp(a.type,'ellipse2D')
                        sel.poly(k).special = EllipseVector2Sym(a.vectors,a.logic);
                    end
                end
            end
        end
    end
    
    % Get
    methods
        function b = vide(sel)
            b = isempty(sel.poly);
        end
        function b = ispoint(sel,tol)
            if ~isscalar(sel.poly)
                error('''ispoint'' method cannot be applied only on non-empty and non-composite selection')
            end
            if nargin<2, tol = 0; end
            switch sel.poly.type
                case {'point1D','point2D','line1D','line2D','poly2D'}
                    points = sel.poly.points;
                    if isempty(points), error programming, end
                    b = all(all(abs(diff(points,1,2))<=tol));
                case 'rect2D'
                    b = all(abs(sel.poly.vectors)<=tol);
                case {'ellipse2D' 'ring2D'}
                    b = all(abs(sel.poly.vectors)<=tol);
                case {'all1D' 'all2D'}
                    b = false;
                otherwise
                    error programming
            end
        end
        function t = get.type(sel)
            if vide(sel)
                t = '';
            else
                t = sel.poly(1).type;
                for i=2:length(sel.poly)
                    if ~strcmp(sel.poly(i).type,t)
                        t = 'mixed';
                        return
                    end
                end
            end
        end
        function m = get.mask(sel)
            if isempty(sel.datasizes)
                m = [];
            else
                m = false([sel.datasizes 1]);
                m(sel.dataind) = true;
            end
        end
    end
    
    % Operations 
    % (note: specialized operations are in functions outside the classdef)
    methods
        function sel2 = copy(sel)
            sel2 = xplr.selectionND('empty',sel.nd);
            sel2.id = sel.id;
            sel2.poly = sel.poly;
            sel2.datasizes = sel.datasizes;
            sel2.dataind = sel.dataind;
        end
        function sel = union(sel1,sel2)
            % function sel = union(sel1,sel2)
            %--
            % BEWARE: this is not symmetric!
            % in particular the ID of the first sel is kept (but apparently, there is
            % even more assymetry?)
            sel = copy(sel1(1));
            seladd = sel1(2:end);
            if nargin==2, seladd = [seladd sel2]; end
            dounion(sel,seladd)
        end
        function dounion(sel1,sel2)
            % function dounion(sel1,sel2)
            %--
            % add the content of sel2 to the content of sel1
            % keep the first ID of sel1, make new second ID
            % sel1 must be a singleton, sel2 does not need to be
            
            % checks
            if any(diff([sel1.nd sel2.nd]))
                error 'dimension mismatch'
            end
            siz = cat(1,sel1.datasizes,sel2.datasizes);
            if any(any(diff(siz)))
                error 'selections don''t have the same data sizes'
            end
            
            % change second part of ID (and count how many unions!)
            ids = cat(1,sel1.id,sel2.id);
            sel1.id(2) = 1 + sum(floor(ids(:,2))) + mod(sum(ids(:,2)),1);
            
            % poly -> specialized functions
            switch sel1.nd
                case 1
                    sel1.poly = union1D([sel1.poly sel2.poly]);
                case 2
                    poly2 = [sel2.poly];
                    for k=1:length(poly2)
                        sel1.poly = union2D(sel1.poly,poly2(k));
                    end
            end
            
            % indices
            sel1.dataind = union(row(sel1.dataind),row(sel2.dataind));
        end
        function sel = substitute(sel1,sel2)
            % function sel = substitute(sel1,sel2)
            %--
            % replace sel by sel2, but keep the first ID and active flag of
            % sel
            
            % multiple object
            if ~isscalar(sel1)
                sel = sel1; % now sel and sel1 contain same objects, but all elements of sel will be changed
                for i=1:length(sel), sel(i) = substitute(sel1(i),sel2(i)); end
                return
            end
            
            % check
            if sel2.nd~=sel1.nd || any(sel2.datasizes~=sel1.datasizes)
                error('dimension mismatch')
            end
            if xor(isequal(sel1.dataind,0),isequal(sel2.dataind,0))
                error('indices computed for only one of the two selections to substitute')
            end
            
            % copy the appropriate parts of sel1 and sel2 (new id will be
            % part from sel1, part from sel2)
            sel = selectionND('empty',sel1.nd);
            sel.id = sel1.id;
            sel.id(2) = sel2.id(2);
            sel.poly = sel2.poly;
            sel.datasizes = sel1.datasizes;
            sel.dataind = sel2.dataind;
        end
        function dosubstitute(sel1,sel2)
            % function sel = substitute(sel1,sel2)
            %--
            % replace content of sel1 by content of sel2
            % keep the first ID of sel1, and the second ID of sel2
            
            % multiple object
            if ~isscalar(sel1)
                for i=1:length(sel1), dosubstitute(sel1(i),sel2(i)); end
                return
            end
            
            % check
            if sel2.nd~=sel1.nd || any(sel2.datasizes~=sel1.datasizes)
                error('dimension mismatch')
            end
            if xor(isequal(sel1.dataind,0),isequal(sel2.dataind,0))
                error('indices computed for only one of the two selections to substitute')
            end
            
            % copy the appropriate parts of sel1 and sel2 (new id will be
            % part from sel1, part from sel2)
            sel1.id(2) = sel2.id(2);
            sel1.poly = sel2.poly;
            sel1.dataind = sel2.dataind; % datasizes is unchanged
        end
        function sel = selaffinity(sel1,mat,varargin)
            % function sel = selaffinity(sel,mat[,datasizes])
            %---
            % mat can be either an affinity matrix, or an affinityND object
            % - in the first case, the ID is unchanged (sel1 and sel
            % represent the same selection, in different referentials);
            % indices are re-computed in sel if the new data sizes are
            % specified as a 3rd argument
            % - in the second case, the second ID number is changed according
            % to the ID of the affinityND object (sel is a transformation
            % of sel1, inside a single same referential); indices are
            % re-computed in sel if they were in sel1
            
            % multiple object
            if ~isscalar(sel1)
                sel = sel1; % now sel and sel1 contain same objects, but elements of sel will be changed
                for i=1:length(sel1), sel(i) = selaffinity(sel1(i),mat,varargin{:}); end 
                return
            end
            
            % note that:
            % - change in ID iff mat is an affinityND object
            % - 'poly' property: fields 'type' and 'logic' remain
            %   unchanged, need to update fields 'points' and 'vectors'
            % - indices must be reset
            sel = copy(sel1);
            
            if isa(mat,'affinityND')
                sel.id(2) = mod(sel.id(2)+mat.id,1);
                mat = mat.mat; % he he... (mat becomes a matrix)
                doindices = true;
                if nargin==3
                    datasizesnew = varargin{1};
                else
                    datasizesnew = sel1.datasizes;
                end
                sel.datasizes = []; % forces new indices computation, see 'ComputeInd' function
            else
                doindices = (nargin==3);
                if doindices, datasizesnew = varargin{1}; end
            end
            
            % check
            if ~all(size(mat)==sel1.nd+1)
                error('wrong size for affinity matrix')
            end
            
            % perform operation 
            for k=1:length(sel.poly)
                a = sel.poly(k);
                xpoints0 = a.points;
                xpoints = [ones(1,size(xpoints0,2)); xpoints0];
                ypoints = mat*xpoints;
                sel.poly(k).points = ypoints(2:end,:);
                switch a.type
                    case {'ellipse2D' 'ring2D'}
                        % this is a bit bad that ellipse main axis vector
                        % and eccentricity cannot be dealt like usual
                        % vectors and logic
                        [sel.poly(k).vectors sel.poly(k).logic sel.poly(k).special] = ...
                            EllipseAffinity(a.vectors,a.logic(1),a.special,mat(2:3,2:3));
                        if strcmp(a.type,'ring2D'), sel.poly(k).logic(2) = a.logic(2); end
                    otherwise
                        sel.poly(k).vectors = mat(2:end,2:end)*a.vectors;
                        % 'logic' remains unchanged
                end
            end
            
            % compute indices
            if doindices && ~isempty(datasizesnew)
                sel = ComputeInd(sel,datasizesnew);
            else
                sel.datasizes = [];
                sel.dataind = 0;
            end
        end
        function sel2 = ComputeInd(sel,datasizes)
            % find indices of an array of size 'datasizes' which are inside
            % the selection described by 'sel'

            % multiple object
            if ~isscalar(sel)
                sel2 = sel;
                for i=1:length(sel), sel2(i) = ComputeInd(sel(i),datasizes); end
                return
            end
            
            % check
            if length(datasizes)~=sel.nd
                error('wrong number of data sizes')
            end
            if isequal(sel.datasizes,datasizes)
                % indices already computed - no need to continue
                sel2 = sel;
                return
            elseif isempty(sel.datasizes)
                % it is ok to have the change to apply on sel as well
                sel2 = sel;
            else
                % the different data sizes for sel might be needed
                % elsewhere: better to make a copy
                sel2 = copy(sel);
            end
            
            % set datasizes
            sel2.datasizes = datasizes;
            
            % specialized functions
            switch sel.nd
                case 1
                    sel2.dataind = indices1D(sel2.poly,datasizes);
                case 2
                    sel2.dataind = indices2D(sel2.poly,datasizes);
            end
        end     
        function sel2 = convert(sel,typenew)
            % function sel2 = convert(sel,'line1D|poly2D')
            %---
            % convert the composite selection into a selection with a
            % unique element
            
            % multiple object
            if ~isscalar(sel)
                sel2 = sel;
                for i=1:length(sel), sel2(i) = convert(sel(i),typenew); end 
                return
            end
            
            % check
            if strcmp(sel.type,typenew), sel2=sel; return, end
            
            % convert
            switch lower(typenew)
                case 'line1d'
                    if sel.nd~=1, error('selection must be 1D'), end
                    polynew = ConvertLine1D(sel.poly);
                case 'poly2d'
                    if sel.nd~=2, error('selection must be 2D'), end
                    polynew = sel.poly; polynew(:) = [];
                    for i=1:length(sel.poly)
                        polynew = union2D(polynew,ConvertPoly2D(sel.poly(i)));
                    end
                otherwise
                    error('conversion to type ''%s'' not implemented',typenew)
            end
            sel2 = copy(sel);
            sel2.poly = polynew;
        end
        function sel2 = convertdim(sel,ndnew,sizesnew)
            % function sel2 = convert(sel,ndnew,sizesnew)
            %---
            % change the dimensionality of a selection; obviously the
            % problem of how to do it is ill-posed
            
            
            % multiple object
            if ~isscalar(sel)
                sel2 = sel;
                if nargin>=3, sizesnew={sizesnew}; else sizesnew={}; end
                for i=1:length(sel), sel2(i) = convert(sel(i),ndnew,sizesnew{:}); end 
                return
            end
            
            % check
            if ndnew==sel.nd, sel2 = sel; return, end
            
            % convert
            switch ndnew
                case 1
                    if sel.nd~=2, error 'case not handled yet', end
                    poly2d = ConvertPoly2D(sel.poly);
                    poly2d = cat(3,poly2d.points);
                    line1d = struct('type','line1D','points',[min(poly2d(1,:)) max(poly2d(1,:))], ...
                        'vectors',[],'logic',[],'special',[]);
                    sel2 = selectionND;
                    sel2.poly = line1d;
                case 2
                    if sel.nd~=1, error 'case not handled yet', end
                    line1d = ConvertLine1D(sel.poly);
                    line1d = cat(3,line1d.points);
                    if nargin>=3, height = sizesnew(2); else height = 1; end
                    rect2d = line1d; rect2d(2,1,:) = .5; rect2d(2,2,:) = height+.5;
                    rect2d = struct('type','rect2D','points',num2cell(rect2d(:,1,:),1), ...
                        'vectors',num2cell(diff(rect2d,1,2),1),'logic',[],'special',[]);
                    sel2 = selectionND;
                    sel2.poly = rect2d;
                otherwise
                    error 'case not handled yet'
            end
            sel2.nd = ndnew;
            sel2.id = sel.id;
            if nargin>=3, sel2 = ComputeInd(sel2,sizesnew); end
        end
    end
    
    % User
    methods
        function [dataind mask] = realworld2dataindices(sel,mat,datasizes)
            % function [dataind mask] = realworld2dataindices(sel,mat,datasizes)
            %---
            % this function returns the indices of data points that are
            % inside a selection expressed in real-world coordinates
            mat = buildMat(mat,1,sel.nd);
            if rank(mat)~=length(mat), error 'ill-defined transformation matrix', end 
            datasel = selaffinity(sel,mat^-1,datasizes);
            dataind = datasel.dataind;
            if nargout>=2
                mask = false([datasizes 1]);
                mask(dataind) = true;
            end
        end
    end
end

%------------
% POLY UNION
%------------

function poly = union1D(poly)

if any(strcmp({poly.type},'all1D'))
    poly = struct('type','all1D','points',[],'vectors',[],'logic',[],'special',[]); 
    return
end
polypoint = poly(strcmp({poly.type},'point1D'));
if ~isempty(polypoint), polypoint = struct('type','point1D','points',unique([polypoint.points]),'vectors',[],'logic',[],'special',[]); end
polyline = poly(strcmp({poly.type},'line1D'));
if ~isempty(polyline), polyline = ConvertLine1D(polyline); end
poly = [polypoint polyline];

end

%---
function poly = union2D(poly,poly2)
% poly2 is a singleton structure

if any(strcmp({poly.type},'all2D'))
    poly = struct('type','all2D','points',[],'vectors',[],'logic',[],'special',[]); 
    return
end
switch poly2.type
    case 'point2D'
        f = find(strcmp({poly.type},'point2D'));
        switch length(f)
            case 0
                poly = [poly poly2];
            case 1
                poly(f).points = [poly(f).points poly2.points];
            otherwise
                error('there should be only one ''point2D'' element in poly')
        end
    case 'line2D'
        error('union of 2D-selection(s) of type ''line2D'' impossible')
    case 'poly2D'
        f = find(strcmp({poly.type},'poly2D'));
        switch length(f)
            case 0
                poly = [poly poly2];
            case 1
                poly(f).points = fn4D_polyunion(poly(f).points,poly2(f).points);
            otherwise
                error('there should be only one ''poly2D'' element in poly')
        end
    case {'rect2D' 'ellipse2D' 'ring2D'}
        poly = [poly poly2];
    otherwise
        error programming
end

end

%---------
% INDICES
%---------

function ind = indices1D(poly,sizes)
ind = zeros(1,0,'uint32');
for k=1:length(poly)
    switch poly(k).type
        case 'point1D'
            ind = [ind round(poly(k).points)]; %#ok<AGROW>
        case 'line1D'
            line = poly(k).points;
            istart = max(1,ceil(line(1)));
            iend   = min(sizes,floor(line(2)));
            ind = [ind istart:iend]; %#ok<AGROW>
        case 'all1D'
            ind = ':';
            return
    end
end
ind(ind<1 | ind>sizes) = [];
ind = unique(ind);
end

%---
function ind = indices2D(poly,sizes)
mask = false([sizes 1]);
for k=1:length(poly)
    points = poly(k).points;
    vects  = poly(k).vectors;
    switch poly(k).type
        case 'point2D'
            np = size(points,2);
            ij = round(points);
            bad = any(ij<1 | ij>repmat(sizes(:),1,np));
            ij(:,bad) = [];
            mask(ij(1,:)+sizes(1)*(ij(2,:)-1)) = true;
        case 'line2D'
            line = points;
            longueur = [0 cumsum(sqrt(sum(diff(line,1,2).^2)))];
            [longueur f] = unique(longueur);
            line = line(:,f);
            longueur2 = [0:.1:longueur(end)-.05 longueur(end)];
            line2 = interp1(longueur,line',longueur2,'pchip')';
            np = size(line2,2);
            ij = round(line2);
            bad = any(ij<1 | ij>repmat(sizes(:),1,np));
            ij(:,bad) = [];
            mask(ij(1,:)+sizes(1)*(ij(2,:)-1)) = true;
        case 'poly2D'
            pp = mypolysplit(points);
            for i=1:length(pp)
                mask = mask | fn_poly2mask(pp{i},sizes);
            end
        case 'rect2D'
            istart = max(1,ceil(points(1)));
            iend   = min(sizes(1),floor(points(1)+vects(1)));
            jstart = max(1,min(ceil(points(2)),round(points(2)+vects(2))));
            jend   = min(sizes(2),max(floor(points(2)+vects(2)),round(points(2))));
            mask(istart:iend,jstart:jend) = true;
        case 'ellipse2D'
            tmp = ConvertPoly2D(poly(k));
            points = tmp.points;
            mask = mask | fn_poly2mask(points,sizes);
        case 'ring2D'
            tmp = ConvertPoly2D(poly(k));
            pp = mypolysplit(tmp.points);
            mask = mask | xor(fn_poly2mask(pp{1},sizes),fn_poly2mask(pp{2},sizes));
        case 'all2D'
            ind = ':';
            return
        otherwise
            error programming
    end
end
ind = uint32(find(mask));
end

%------------
% CONVERSION
%------------

function poly2 = ConvertLine1D(poly)
    lines = zeros(0,2);
    for k=1:length(poly)
        switch poly(k).type
            case 'line1D'
                lines(end+1,:) = poly(k).points; %#ok<AGROW>
            case 'point1D'
                pointsk = poly(k).points;
                linesk = fn_add(round(pointsk)',[-.5 .5]);
                lines = [lines; linesk]; %#ok<AGROW>
            case 'all1D'
                disp 'converting ''all1D'' selection by a [-1e30 1e30] line'
                lines = [-1 1]*1e30;
                break
            otherwise
                error programming
        end
    end
    merge = xor( bsxfun(@lt,lines(:,1),lines(:,2)'), bsxfun(@lt,lines(:,2),lines(:,1)') );
    for i=1:size(lines,1), merge(i,i) = false; end
    while any(merge(:))
        [i j] = find(merge,1,'first');
        lines(i,:) = [min(lines([i j],1)) max(lines([i j],2))];
        merge(i,:) = merge(i,:) | merge(j,:);
        merge(:,i) = merge(:,i) | merge(:,j);
        merge(i,i) = false;
        lines(j,:) = [];
        merge(j,:) = [];
        merge(:,j) = [];
    end
    lines = num2cell(lines,2)';
    poly2 = struct('type','line1D','points',lines,'vectors',[],'logic',[],'special',[]);
end

function poly2 = ConvertPoly2D(poly)
    if ~isscalar(poly), error('poly must be a scalar structure here'), end
    switch poly.type
        case 'poly2D'
            poly2 = poly;
            return
        case 'point2D'
            points = round(poly.points);
            points2 = zeros(2,0);
            for i=1:size(points,2)
                points2 = fn4D_polyunion(points2, ...
                    fn_add(points(:,i),[-.5 -.5 .5 .5 -.5; -.5 .5 .5 -.5 -.5]));
            end
        case 'rect2D'
            points2 = fn_add(poly.points,fn_mult(poly.vectors,[0 0 1 1 0; 0 1 1 0 0]));
        case {'ellipse2D' 'ring2D'}
            c = poly.points;
            u = poly.vectors;
            e = poly.logic(1);
            phi = linspace(0,2*pi,20);
            udata = cos(phi);
            vdata = e*sin(phi);
            if strcmp(poly.type,'ring2D')
                r = poly.logic(2);
                udata = [udata NaN r*udata];
                vdata = [vdata NaN r*vdata];
            end
            points2 = fn_add(c,fn_mult(u,udata)+fn_mult([u(2);-u(1)],vdata));
        case 'line2D'
            points2 = poly.points;
        case 'all2D'
            disp 'converting ''all2D'' selection by a [-1e30 1e30] square'
            points2 = [-1 -1 1 1 -1; -1 1 1 -1 -1]*1e30;
        otherwise
            error programming
    end
    poly2 = struct('type','poly2D','points',points2,'vectors',zeros(2,0), ...
        'logic',[],'special',[]);
end

%------------
% ELLIPSE
%------------

% Ellipse defined either by center, main radius vector and eccentricity, or
% by its bilinear equation (x-c)'A(x-c) = 1,

function [u, e, A] = EllipseAffinity(u,e,A,M)

% ellipse equation becomes, for y=Mx: (y-Mc)'(M^-1' A M^-1)(y-Mc) = 1
M1 = M^-1;
A = M1'*A*M1;
[u, e] = EllipseSym2Vector(A);

end

%---
function A = EllipseVector2Sym(u,e)

% (U,c) is the referential of the ellipse, x->y=U'(x-c) returns coordinates
% in this referential, in which the ellipse equation is 
% y(1)^2 + y(2)^2/e^2 = r^2
r = norm(u);
u = u/r;
U = [[-u(2); u(1)] u];        
A = U*(diag([1/e^2 1])/r^2)*U';

end

%---
function [u e] = EllipseSym2Vector(A)

% eigenvalue decomposition returns U, r and e as above
[U D] = svd(A); % better use svd than eig because output is real even if A is not exactly symmetric
r = D(2,2)^-(1/2);
e = D(1,1)^-(1/2) / r;
u = r*U(:,2);

end

%---
function pp = mypolysplit(points)

ksep = find(any(isnan(points),1));
n = length(ksep)+1;
if n==1
    pp = {points};
else
    pp = cell(1,n); ok = false(1,n);
    ksep = [0 ksep size(points,2)+1];
    for i=1:n
        pp{i} = points(:,ksep(i)+1:ksep(i+1)-1);
        ok(i) = ~isempty(pp{i});
    end
    pp = pp(ok);
end

end
