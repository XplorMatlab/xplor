classdef selectionnd < xplr.graphnode
    % function sel = selectionnd('type',data[,sizes])
    %---
    % selectionnd class defines selection in an formal manner, i.e. by
    % defining the shape that makes the selection. This shape can
    % thereafter be converted to the indices of the data points they cover
    % once the sizes of the data are provided. Yet shape do not rely on
    % specific data sizes, and geometric operations can be applied to them
    % (affine transformation, union, etc.).
    %     
    % Available shape type are:
    % - 1D shapes: 'point1D', 'line1D'
    % - 2D shapes: 'point2D', 'line2D', 'poly2D', 'rect2D', 
    %              'ellipse2D', 'ring2D'
    % - 'indices': those are not shapes but directly data indices,
    %   therefore data sizes need to be provided and no affine
    %   transformation can be applied on them
    % - ND shapes: 'product' combines selection in multiple dimensions,
    %              'emptyND', 'allND' (replace 'N' by the number of
    %              dimensions)
    %
    % If 'sizes' is specified, indices are computed according to this data
    % sizes. 'sizes' is mandatory for 'indices' selections.
    % 
    % 'selectionnd' is a handle class to make a = b faster (no copy), but
    % be careful when using it.
    % 
    % See xplr.shelectionshape for details on the 'data' argument.
    % 
    % See also xplr.selectionshape

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
        shapes = xplr.selectionshape.empty(1,0);
        datasizes = [];
        dataind = []; % to not be mistaken with zeros(1,0)
    end
    properties (Dependent, SetAccess='private')
        type % either a possible value for shapes.type (e.g. point1D, ellipse2D), or 'mixed' if several types are present
        mask
    end
    
    % Constructor + Load + Display
    methods
        function sel = selectionnd(type,data,sizes)
            % function sel = selectionnd('type',data[,sizes])
            %---
            % selection ID is made of 2 numbers: the first remains the same as long as
            % the selection is living; the second is changed each time the seleciton is
            % modified
            
            % empty selection (invalid)
            if nargin==0
                sel.nd = 0;
                return
            end
            
            % random ID
            sel.id = rand(1,2);
            
            % number of dimension
            if fn_ismemberstr(type, {'empty' 'all'})
                % the syntax selectionnd('empty',nd) is not public but
                % corresponds to the internal encoding of type
                sel.nd = data;
                data = [];
            elseif regexp(type, '^(empty|all)\d+D$')
                [type, nd] = fn_regexptokens(type, '^(empty|all)(\d+)D$');
                sel.nd = str2double(nd);
            elseif contains(type,'1D')
                sel.nd = 1;
            elseif contains(type,'2D')
                sel.nd = 2;
            elseif strcmp(type, 'indices')
                if nargin<=3
                    error 'Data sizes must be provided for ''indices'' selectionnd type'
                end
                sel.nd = length(sizes);
            else
                error('unknown selection type ''%s''',type)
            end
            
            % active by default
            sel.active = true;
            
            % shape
            if ~strcmp(type,'empty')
                sel.shapes = xplr.selectionshape(type,data);
            end
                        
            % indices
            if nargin==3
                sel = ComputeInd(sel,sizes);
            end
        end
        function disp(sel)
            warning('off','MATLAB:structOnObject')
            if isempty(sel)
                disp('empty selectionnd object')
            elseif ~isscalar(sel)
                s = size(sel);
                str = cellstr(num2str(s'))';
                [str{2,:}] = deal('x');
                str = [str{1:end-1} ' selectionnd object'];
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
            % previous version might not have a shapes.special field
            if ~isfield(sel.shapes,'special')
                if isempty(sel.shapes)
                    error('programming: a selection cannot have an empty shapes')
                end
                sel.shapes(1).special = [];
                for k=1:length(sel.shapes)
                    a = sel.shapes(k);
                    if strcmp(a.type,'ellipse2D')
                        sel.shapes(k).special = EllipseVector2Sym(a.vectors,a.logic);
                    end
                end
            end
        end
    end
    
    % Misc
    methods
        function b = vide(sel)
            b = isempty(sel.shapes);
        end
        function b = ispoint(sel,tol)
            if ~isscalar(sel.shapes)
                error('''ispoint'' method cannot be applied only on non-empty and non-composite selection')
            end
            if nargin<2, tol = 0; end
            b = ispoint(sel.shapes,tol);
        end
        function t = get.type(sel)
            if vide(sel)
                t = '';
            else
                t = sel.shapes(1).type;
                for i=2:length(sel.shapes)
                    if ~strcmp(sel.shapes(i).type,t)
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
            sel2 = xplr.selectionnd('empty',sel.nd);
            sel2.id = sel.id;
            sel2.shapes = sel.shapes;
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
            
            % shapes
            sel1.shapes = simplify([sel1.shapes sel2.shapes]);
            
            % indices
            if ~isempty(siz)
                if all(strcmp(sel1.shapes.type,'indices'))
                    sel1.dataind = unique([sel1.shapes.special]);
                else
                    sel1.dataind = union(sel1.dataind,sel2.dataind);
                end
            end
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
            if xor(isequal(sel1.dataind,[]),isequal(sel2.dataind,[]))
                error('indices computed for only one of the two selections to substitute')
            end
            
            % copy the appropriate parts of sel1 and sel2 (new id will be
            % part from sel1, part from sel2)
            sel = selectionnd('empty',sel1.nd);
            sel.id = sel1.id;
            sel.id(2) = sel2.id(2);
            sel.shapes = sel2.shapes;
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
            if xor(isequal(sel1.dataind,[]),isequal(sel2.dataind,[]))
                error('indices computed for only one of the two selections to substitute')
            end
            
            % copy the appropriate parts of sel1 and sel2 (new id will be
            % part from sel1, part from sel2)
            sel1.id(2) = sel2.id(2);
            sel1.shapes = sel2.shapes;
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
            % - 'shapes' property: fields 'type' and 'logic' remain
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
            sel.shapes = affinity(sel.shapes,mat);
            
            % compute indices
            if doindices && ~isempty(datasizesnew)
                sel = ComputeInd(sel,datasizesnew);
            else
                sel.datasizes = [];
                sel.dataind = [];
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
            
            % indices
            if all(strcmp(sel2.type),'indices')
                sel2.dataind = unique([sel2.shapes.special]);
            else
                switch sel.nd
                    case 1
                        sel2.dataind = indices1D(sel2.shapes,datasizes);
                    case 2
                        sel2.dataind = indices2D(sel2.shapes,datasizes);
                end
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
                    shapesnew = ConvertLine1D(sel.shapes);
                case 'poly2d'
                    if sel.nd~=2, error('selection must be 2D'), end
                    shapesnew = xplr.selectionshape.empty(1,0);
                    for i=1:length(sel.shapes)
                        shapesnew = union2D(shapesnew,ConvertPoly2D(sel.shapes(i)));
                    end
                otherwise
                    error('conversion to type ''%s'' not implemented',typenew)
            end
            sel2 = copy(sel);
            sel2.shapes = shapesnew;
        end
    end
    
    % User
    methods
        function [dataind, mask] = realworld2dataindices(sel,mat,datasizes)
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

