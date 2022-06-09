classdef SelectionND < xplr.Object
    % function sel = SelectionND('type', data[, sizes])
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
    % - 'shapeND' (replace 'N' by number of dimensions): data must already
    %   be an xplr.SelectionShape object 
    %
    % If 'sizes' is specified, indices are computed according to this data
    % sizes. 'sizes' is mandatory for 'indices' selections.
    % 
    % 'SelectionND' is a handle class to make a = b faster (no copy), but
    % be careful when using it.
    % 
    % See xplr.SelectionShape for details on the 'data' argument.
    % 
    % See also xplr.SelectionShape

    properties (SetAccess='private')
        nd
        shapes = xplr.SelectionShape.empty(1,0);
        data_sizes = [];
    end
    properties (SetAccess='private', Transient)
        data_ind = [];  % not computed yet, to not be mistaken with zeros(1,0) = empty selection
        polygon = []; % lines(1D)/polygon(2D) to display, not computed yet, to not be mistaken with zeros(2,0) = empty selection
    end
    properties (Dependent, SetAccess='private', Transient)
        type % either a possible value for shapes.type (e.g. point1D, ellipse2D), or 'mixed' if several types are present
        mask
    end
    
    % Constructor + Load + Display
    methods
        function sel = SelectionND(type, data, sizes)
            % function sel = SelectionND('type', data[, sizes])
            
            % initialize a single object or an array thereof, invalid(s)
            % for the moment
            if nargin==0
                sel.nd = 0;
                return
            elseif isnumeric(type)
                n = type;
                if n == 0
                    sel = xplr.SelectionND.empty(1,0);
                else
                    sel(n) = xplr.SelectionND;
                end
                return
            end
            
            % number of dimension
            if nargin == 3, sizes = brick.row(sizes); end
            if brick.ismemberstr(type, {'empty', 'all'})
                % the syntax SelectionND('empty', nd) is not public but
                % corresponds to the internal encoding of type
                sel.nd = data;
                data = [];
            elseif regexpi(type, '^(empty|all|shape)\d+D$')
                [type, nd] = brick.regexptokens(type, '^(empty|all|shape)(\d+)D$');
                type = lower(type);
                nd = upper(nd);
                sel.nd = str2double(nd);
            elseif strfind(upper(type), '1D')
                sel.nd = 1;
            elseif strfind(upper(type), '2D')
                sel.nd = 2;
            elseif strcmp(type, 'indices')
                % two syntaxes are accepted:
                % xplr.SelectionND('indices', indices,datasize)
                % xplr.SelectionND('indices', {datasize indices})
                if nargin==2 && iscell(data)
                    sizes = brick.row(data{1});
                elseif nargin==3 && ~iscell(data)
                    indices = brick.row(data);
                    data = {sizes, indices};
                else
                    error argument
                end
                sel.nd = length(sizes);
            else
                error('unknown selection type ''%s''', type)
            end
            
            % shape
            switch type
                case 'empty'
                    % no shape
                case 'shape'
                    sel.shapes = data;
                otherwise
                    if isempty(data)
                        % happens sometimes: points selection with no
                        % poins, empty region, etc.
                        type = 'empty';
                    else
                        sel.shapes = xplr.SelectionShape(type, data);
                    end
            end
                        
            % indices
            if nargin==3
                sel = compute_indices(sel, sizes);
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
                brick.strucdisp(struct(sel))
            end
            warning('on', 'MATLAB:structOnObject')
        end
    end
    
    % Info
    methods
        function b = vide(sel)
            b = isempty(sel.shapes);
        end
        function b = is_point(sel, tol)
            if ~isscalar(sel.shapes)
                error('''is_point'' method cannot be applied only on non-empty and non-composite selection')
            end
            if nargin<2, tol = 0; end
            b = is_point(sel.shapes,tol);
        end
        function check_point(sel, tol)
            if sel.is_point(tol)
                sel.shapes = sel.shapes.to_point();
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
            sel2 = xplr.SelectionND('empty', sel.nd);
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
            do_union(sel, seladd)
        end
        function do_union(sel1, sel2)
            % function do_union(sel1, sel2)
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
        function sel2 = apply_affinity(sel1, affinity, data_sizesnew)
            % function sel2 = apply_affinity(sel, mat[, data_sizesnew])
            %---
            % affinity is an xplr.AffinityND object
            
            if isempty(sel1)
                sel2 = xplr.SelectionND.empty(size(sel1));
                return
            elseif nargin<3
                data_sizesnew = sel1.data_sizes;
            end
            data_sizesnew = brick.row(data_sizesnew);
            
            % multiple object
            if ~isscalar(sel1)
                sel2 = sel1; % now sel and sel1 contain same objects, but elements of sel will be changed
                for i=1:length(sel1), sel2(i) = apply_affinity(sel1(i), affinity, data_sizesnew); end
                return
            end
            
            sel2 = xplr.SelectionND('empty', sel1.nd);
            
            % perform operation 
            if ~isa(affinity, 'xplr.AffinityND')
                error('argument ''afinity'' is expected to be an affinitynd instance')
            end
            sel2.shapes = apply_affinity(sel1.shapes, affinity);
            
            % compute indices
            if ~isempty(data_sizesnew)
                sel2.data_sizes = []; % forces new indices computation even if data_sizes did not change, see 'compute_indices' function
                sel2 = compute_indices(sel2, data_sizesnew);
            else
                sel2.data_sizes = [];
                sel2.data_ind = [];
            end
            
            % invalidate polygon
            sel2.polygon = [];
            
        end
        function sel2 = compute_indices(sel, data_sizes)
            % find indices of an array of size 'data_sizes' which are inside
            % the selection described by 'sel'

            % multiple object
            if ~isscalar(sel)
                sel2 = sel;
                for i=1:length(sel), sel2(i) = compute_indices(sel(i), data_sizes); end
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
                        sel2.data_ind = sel2.shapes.indices_1D(data_sizes);
                    case 2
                        sel2.data_ind = sel2.shapes.indices_2D(data_sizes);
                end
            end
        end     
        function polygon = get.polygon(sel)
            if isempty(sel.polygon)
                switch sel.nd
                    case 1
                        sel.polygon = sel.shapes.convert_line_1D();
                    case 2
                        sel.polygon = sel.shapes.convert_poly_2D();
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
                sel2 = xplr.SelectionND(length(sel));
                for i = 1:length(sel)
                    sel2(i) = convert(sel(i), type, data_sizes);
                end
                return
            end
            
            % conversion
            switch type
                case 'indices'
                    if nargin<3, data_sizes = sel.data_sizes; end
                    sel2 = compute_indices(sel, data_sizes);
                    if ~strcmp(sel2.type, 'indices')
                        sel2 = xplr.SelectionND('indices', {data_sizes, sel.data_ind});
                    end
                otherwise
                    error 'invalid correction type'
            end
        end
    end
end
