classdef AffinityND < handle
    % function mov = AffinityND(linear_part, offset)
    % function mov = AffinityND('translate1D|2D', offset)
    % function mov = AffinityND('scale1D|2D', scale)
   
    properties (SetAccess='private')
        nd
    end
    properties (Access='private')
        mat
    end
    properties (Dependent)
        linear_part
        offset
    end
    
    % Constructor + Load + Display
    methods
        function mov = AffinityND(varargin)
            if isnumeric(varargin{1}) || isduration(varargin{1})
                [linear_part, offset] = deal(varargin{:});
                mov.nd = length(offset);
                if mov.nd > 1 && isvector(linear_part)
                    linear_part = diag(linear_part);
                end
                mov.mat = [1, zeros(1,mov.nd); offset(:), linear_part];
            else
                [type, data] = deal(varargin{:});
                % number of dimensions
                switch type
                    case {'translate1D', 'scale1D'}
                        mov.nd = 1;
                    case {'translate2D', 'scale2D'}
                        mov.nd = 2;
                    otherwise
                        error('unknown affinity type ''%s''', type)
                end
                
                % affinity matrix
                switch type
                    case 'translate1D'
                        if ~isscalar(data), error('wrong translation data'), end
                        mov.mat = [1, 0; data, 1];
                    case 'scale1D'
                        if ~isscalar(data), error('wrong scaling factor'), end
                        mov.mat = diag([1, data]);
                    case 'translate2D'
                        if ~isvector(data) || length(data)~=2, error('wrong translation data'), end
                        mov.mat = [1, 0, 0; data(:), eye(2)];
                    case 'scale2D'
                        if isscalar(data), data = [data, data]; end
                        if ~isvector(data) || length(data)~=2, error('wrong scaling data'), end
                        mov.mat = diag([1, data(:)']);
                end
            end
        end
        function x = get.linear_part(mov)
            x = mov.mat(2:end, 2:end);
        end
        function x = get.offset(mov)
            x = mov.mat(2:end, 1);
        end
        function disp(impossible_name__) %#ok<MANU>
            warning('off', 'MATLAB:structOnObject')
            varname = inputname(1);
            eval([varname, ' = struct(impossible_name__);'])
            brick.structdisp(varname)
            warning('on', 'MATLAB:structOnObject')
        end
    end
    
    % Operations 
    methods
        function points = move_points(mov, points)
            points = brick.add(mov.mat(2:end, 1), mov.mat(2:end, 2:end)*points);
        end
        function vectors = move_vectors(mov, vectors)
            vectors = mov.mat(2:end, 2:end)*vectors;
        end
        function aff = move_affinity(mov, aff)
            % function mov = movaffinity(mov1, mat)
            %---
            % scheme prevents headache: 
            %
            %   shape1, ref1   -- mov --> shape1, ref2
            %        |                        |
            %       aff                      aff (new)
            %        |                        |
            %        V                        V
            %   shape2, ref1  -- mov --> shape2, ref2
            if isa(aff, 'xplr.AffinityND')
                aff.mat = mov.mat*aff.mat*mov.mat^-1;
            else
                aff = mov.mat*aff*mov.mat^-1;
            end
        end
        function selection = move_selection(mov, selection)
            if ~isa(selection, 'xplr.SelectionND'), error 'input ''selection'' must be an xplr.SelectionND object', end
            selection = selection.apply_affinity(mov);
        end
    end
end

