classdef affinitynd < handle
    % function mov = affinitynd('type',data)
    %---
    % Type can be 'translate2D', 'scale2D'.
   
    properties (SetAccess='private')
        nd
        mat
    end
    
    % Constructor + Load + Display
    methods
        function mov = affinitynd(type,data)
            % number of dimensions
            switch type
                case {'translate1D','scale1D'}
                    mov.nd = 1;
                case {'translate2D','scale2D'}
                    mov.nd = 2;
                otherwise
                    error('unknown affinity type ''%s''',type)
            end
            
            % affinity matrix
            switch type
                case 'translate1D'
                    if ~isscalar(data), error('wrong translation data'), end
                    mov.mat = [1 0 ; data 1];
                case 'scale1D'
                    if ~isscalar(data), error('wrong scaling factor'), end
                    mov.mat = diag([1 data]);
                case 'translate2D'
                    if ~isvector(data) || length(data)~=2, error('wrong translation data'), end
                    mov.mat = [1 0 0; data(:) eye(2)];
                case 'scale2D'
                    if isscalar(data), data = [data data]; end
                    if ~isvector(data) || length(data)~=2, error('wrong scaling data'), end
                    mov.mat = diag([1 data(:)']);
            end
        end
        function disp(impossible_name__) %#ok<MANU>
            warning('off','MATLAB:structOnObject')
            varname = inputname(1);
            eval([varname ' = struct(impossible_name__);'])
            fn_structdisp(varname)
            warning('on','MATLAB:structOnObject')
        end
    end
    
    % Operations 
    methods
        function points = move_points(mov, points)
            points = fn_add(mov.mat(2:end,1), mov.mat(2:end,2:end)*points);
        end
        function vectors = move_vectors(mov, vectors)
            vectors = mov.mat(2:end,2:end)*vectors;
        end
        function aff = move_affinity(mov, aff)
            % function mov = movaffinity(mov1,mat)
            %---
            % scheme prevents headache: 
            %
            %   shape1,ref1   -- mov --> shape1,ref2
            %        |                        |
            %       aff                      aff (new)
            %        |                        |
            %        V                        V
            %   shape2, ref1  -- mov --> shape2,ref2
            if isa(aff,'xplr.affinitynd')
                aff.mat = mov.mat*aff.mat*mov.mat^-1;
            else
                aff = mov.mat*aff*mov.mat^-1;
            end
        end
    end
end

