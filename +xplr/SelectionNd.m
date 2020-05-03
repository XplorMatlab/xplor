classdef SelectionNd < xplr.GraphNode
    % function sel = SelectionNd('type', data[, sizes])
    %---
    % selection_nd class defines selection in an formal manner, i.e. by
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
    % 'selection_nd' is a handle class to make a = b faster (no copy), but
    % be careful when using it.
    % 
    % See xplr.SelectionShape for details on the 'data' argument.
    % 
    % See also xplr.SelectionShape

    properties (SetAccess='private')
        nd
        shapes = xplr.SelectionShape.empty(1,0);
        data_sizes = [];
        data_ind = [];  % not computed yet, to not be mistaken with zeros(1,0) = empty selection
        polygon = []; % lines(1D)/polygon(2D) to display, not computed yet, to not be mistaken with zeros(2,0) = empty selection
    end
    properties (Dependent, SetAccess='private')
        type % either a possible value for shapes.type (e.g. point1D, ellipse2D), or 'mixed' if several types are present
        mask
    end
    
    % Constructor + Load + Display
    methods
        function sel = SelectionNd(type, data, sizes)
            % function sel = SelectionNd('type', data[, sizes])
            
            % initialize a single object or an array thereof, invalid(s)
            % for the moment
            if nargin==0
                sel.nd = 0;
                return
            elseif isnumeric(type)
                n = type;
                if n == 0
                    sel = xplr.SelectionNd.empty(1,0);
                else
                    sel(n) = xplr.SelectionNd;
                end
                return
            end
            
            % number of dimension
            if fn_ismemberstr(type, {'empty', 'all'})
                % the syntax SelectionNd('empty', nd) is not public but
                % corresponds to the internal encoding of type
                sel.nd = data;
                data = [];
            elseif regexp(type, '^(empty|all)\d+D$')
                [type, nd] = fn_regexptokens(type, '^(empty|all)(\d+)D$');
                sel.nd = str2double(nd);
            elseif strfind(type, '1D')
                sel.nd = 1;
            elseif strfind(type, '2D')
                sel.nd = 2;
            elseif strcmp(type, 'indices')
                % two syntaxes are accepted:
                % xplr.SelectionNd('indices', indices,datasize)
                % xplr.SelectionNd('indices', {datasize indices})
                if nargin==2 && iscell(data)
                    sizes = data{1};
                elseif nargin==3 && ~iscell(data)
                    indices = row(data);
                    data = {sizes, indices};
                else
                    error argument
                end
                sel.nd = length(sizes);
            else
                error('unknown selection type ''%s''', type)
            end
            
            % shape
            if ~strcmp(type, 'empty')
                sel.shapes = xplr.SelectionShape(type, data);
            end
                        
            % indices
            if nargin==3
                sel = ComputeInd(sel, sizes);
            end
        end
        function disp(sel)
            warning('off', 'MATLAB:structOnObject')
            if isempty(sel)
                disp('empty selection_nd object')
            elseif ~isscalar(sel)
                s = size(sel);
                str = cellstr(num2str(s'))';
                [str{2, :}] = deal('x');
                str = [str{1:end-1}, ' selection_nd object'];
                fprintf('%s\n\nProperties:\n', str)
                F = fieldnames(sel);
                nF = length(F);
                for k=1:nF
                    fprintf('    %s\n', F{k});
                end
            else
                strucdisp(struct(sel))
            end
            warning('on', 'MATLAB:structOnObject')
        end
    end
    
    % Info
    methods
        function b = vide(sel)
            b = isempty(sel.shapes);
        end
        function b = ispoint(sel, tol)
            if ~isscalar(sel.shapes)
                error('''ispoint'' method cannot be applied only on non-empty and non-composite selection')
            end
            if nargin<2, tol = 0; end
            b = ispoint(sel.shapes,tol);
        end
        function checkpoint(sel, tol)
            if sel.ispoint(tol)
                sel.shapes = sel.shapes.topoint();
            end
        end
        function t = get.type(sel)
            if vide(sel)
                t = '';
            else
                t = sel.shapes(1).type;
                for i=2:length(sel.shapes)
                    if ~strcmp(sel.shapes(i).type, t)
                        t = 'mixed';
                        return
                    end
                end
            end
        end
        function m = get.mask(sel)
            if isempty(sel.data_sizes)
                m = [];
            else
                m = false([sel.data_sizes 1]);
                m(sel.data_ind) = true;
            end
        end
    end
    
    % Operations
    % (note: specialized operations are in functions outside the classdef)
    methods
        function sel2 = copy(sel)
            sel2 = xplr.SelectionNd('empty', sel.nd);
            sel2.shapes = sel.shapes;
            sel2.data_sizes = sel.data_sizes;
            sel2.data_ind = sel.data_ind;
            sel2.polygon = sel.polygon;
        end
        function sel = union(sel1, sel2)
            % function sel = union(sel1,sel2)
            sel = copy(sel1(1));
            seladd = sel1(2:end);
            if nargin == 2, seladd = [seladd, sel2]; end
            dounion(sel, seladd)
        end
        function dounion(sel1, sel2)
            % function dounion(sel1, sel2)
            %--
            % add the content of sel2 to the content of sel1
            
            % checks
            if any(diff([sel1.nd, sel2.nd]))
                error 'dimension mismatch'
            end
            siz = cat(1, sel1.data_sizes, sel2.data_sizes);
            if any(any(diff(siz)))
                error 'selections don''t have the same data sizes'
            end
            
            % shapes
            sel1.shapes = simplify([sel1.shapes, sel2.shapes]);
            
            % indices
            if ~isempty(siz)
                if all(strcmp(sel1.shapes.type, 'indices'))
                    sel1.data_ind = unique([sel1.shapes.special]);
                else
                    sel1.data_ind = union(sel1.data_ind, sel2.data_ind);
                end
            end
            
            % invalidate polygon
            sel1.polygon = [];
        end
        function sel2 = applyaffinity(sel1, affinity, data_sizesnew)
            % function sel2 = applyaffinity(sel, mat[, data_sizesnew])
            %---
            % affinity is an xplr.AffinityNd object
            
            if isempty(sel1)
                sel2 = xplr.SelectionNd.empty(size(sel1));
                return
            elseif nargin<3
                data_sizesnew = sel1.data_sizes;
            end
            
            % multiple object
            if ~isscalar(sel1)
                sel2 = sel1; % now sel and sel1 contain same objects, but elements of sel will be changed
                for i=1:length(sel1), sel2(i) = applyaffinity(sel1(i), affinity, data_sizesnew); end
                return
            end
            
            sel2 = xplr.SelectionNd('empty', sel1.nd);
            
            % perform operation 
            if ~isa(affinity, 'xplr.AffinityNd')
                error('argument ''afinity'' is expected to be an affinitynd instance')
            end
            sel2.shapes = applyaffinity(sel1.shapes, affinity);
            
            % compute indices
            if ~isempty(data_sizesnew)
                sel2.data_sizes = []; % forces new indices computation even if data_sizes did not change, see 'ComputeInd' function
                sel2 = ComputeInd(sel2, data_sizesnew);
            else
                sel2.data_sizes = [];
                sel2.data_ind = [];
            end
            
            % invalidate polygon
            sel2.polygon = [];
            
        end
        function sel2 = ComputeInd(sel, data_sizes)
            % find indices of an array of size 'data_sizes' which are inside
            % the selection described by 'sel'

            % multiple object
            if ~isscalar(sel)
                sel2 = sel;
                for i=1:length(sel), sel2(i) = ComputeInd(sel(i), data_sizes); end
                return
            end
            
            % check
            if length(data_sizes) ~= sel.nd
                error('wrong number of data sizes')
            end
            if isequal(sel.data_sizes, data_sizes)
                % indices already computed - no need to continue
                sel2 = sel;
                return
            elseif isempty(sel.data_sizes)
                % it is ok to have the change to apply on sel as well
                sel2 = sel;
            else
                % the different data sizes for sel might be needed
                % elsewhere: better to make a copy
                sel2 = copy(sel);
            end
            
            % set data_sizes
            sel2.data_sizes = data_sizes;
            
            % indices
            if all(strcmp(sel2.type, 'indices'))
                data = cat(1, sel2.shapes.special); % first column: sizes, second column: indices
                sel2.data_ind = unique([data{:, 2}]);
            else
                switch sel.nd
                    case 1
                        sel2.data_ind = indices1D(sel2.shapes, data_sizes);
                    case 2
                        sel2.data_ind = indices2D(sel2.shapes, data_sizes);
                end
            end
        end     
        function polygon = get.polygon(sel)
            if 1 %isempty(sel.polygon)
                switch sel.nd
                    case 1
                        sel.polygon = ConvertLine1D(sel.shapes);
                    case 2
                        sel.polygon = ConvertPoly2D(sel.shapes);
                    otherwise
                        error 'not implemented yet'
                end
            end
            polygon = sel.polygon;
        end
        function sel2 = convert(sel, type, data_sizes)
            % function sel2 = convert(sel, 'indices', data_sizes)
            
            % multiple selections
            if ~isscalar(sel)
                sel2 = xplr.SelectionNd(length(sel));
                for i = 1:length(sel)
                    sel2(i) = convert(sel(i), type, data_sizes);
                end
                return
            end
            
            % conversion
            switch type
                case 'indices'
                    if nargin<3, data_sizes = sel.data_sizes; end
                    sel2 = ComputeInd(sel, data_sizes);
                    if ~strcmp(sel2.type, 'indices')
                        sel2 = xplr.SelectionNd('indices', {data_sizes, sel.data_ind});
                    end
                otherwise
                    error 'invalid correction type'
            end
        end
    end
end
