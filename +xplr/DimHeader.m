classdef DimHeader < xplr.Header
    % function H = dim_header(header arguments...)
    % function H = dim_header(header object[, dim_id])
    %---
    % The dim_header class inherits from HEADER and adds the 'dim_id'
    % property, which is a random number serving to uniquely identify
    % different dimensions inside an xdata object, even if their headers
    % are rigorously the same.
    %
    % See also xplr.header
    
    properties
        dim_id
    end
    
    methods
        function H = DimHeader(varargin)
            % header part
            do_copy = (nargin >= 1 && isa(varargin{1}, 'xplr.Header'));
            if do_copy
                arg = {};
            else
                % construct header with xplr.header syntax
                arg = varargin;
            end
            H = H@xplr.Header(arg{:});
            if do_copy
                % copy header
                H1 = varargin{1};
                if ~isscalar(H1)
                    % pre-allocate
                    H(numel(H1)) = xplr.DimHeader();
                    H = reshape(H, size(H1));
                end
                H.copy_in(H1);
            end
            
            % generate a unique dimension ID
            if do_copy && nargin >= 2
                [H.dim_id] = brick.dealc(varargin{2});
            else
                H.changedim_id()
            end
        end
        function changedim_id(H)
            % generate a new, unique dimension ID; this method will be
            % called in particular when a header will be duplicated (and
            % potentially modified) for a new usage
            for i = 1:length(H)
                H(i).dim_id = rand;
            end
        end
        function disp(H)
            disp@xplr.Header(H)
            fprintf('\b    dim_id: %.4f\n\n', H.dim_id);
        end
        function [dim, dim_id] = dimension_number_and_id(H, d)
            % function [dim, dim_id] = dimension_number_and_id(H,d)
            %---
            % Convert any of dimension numbers, identifiers or labels
            % to both dimension numbers and identifiers.
            % Returns 0 for unfound dimensions.
            
            % Special cases
            if isempty(d)
                if iscell(d)
                    [dim, dim_id] = deal(cell(1,0));
                else
                    [dim, dim_id] = deal(zeros(1,0));
                end
                return
            elseif ischar(d) || (iscell(d) && ischar(d{1}))
                % First convert labels to dimension numbers
                if ~iscell(d), d = {d}; end
                n = length(d);
                dim = zeros(1, n);
                labels = {H.label};
                for i = 1:n
                    dim_i = brick.find(d{i}, labels, 'first');
                    if ~isempty(dim_i)
                        dim(i) = dim_i;
                    end
                end
                d = dim;
            elseif iscell(d)
                % Multiple outputs
                n = length(d);
                [dim, dim_id] = deal(cell(1, n));
                for i = 1:n
                    [dim{i}, dim_id{i}] = H.dimension_number_and_id(d{i}); 
                end
                return
            end
            
            % Convert between dimension numbers and identifiers
            if d(1) < 1
                % identifier -> number
                dim_id = d;
                n = length(dim_id);
                dim = zeros(1, n);
                for i =1:n
                    dim_i = find([H.dim_id] == dim_id(i), 1, 'first');
                    if ~isempty(dim_i)
                        dim(i) = dim_i;
                    end
                end
            else
                % number -> identifier
                dim = d;
                dim_id = [H(dim).dim_id];
                if length(dim_id) < length(d)
                    [dim, dim_id] = deal([]);
                    return
                end
            end
        end
        function dim_id = dimension_id(H, d)
            [~, dim_id] = H.dimension_number_and_id(d);
        end
        function dim = dimension_number(H, d)
            [dim, ~] = H.dimension_number_and_id(d);
        end
        function label = dimension_label(H, d)
            [dim, ~] = H.dimension_number_and_id(d);
            if isempty(dim)
                label = [];
            else
                label = H(dim).label;
            end
        end
        function head = headerByID(H, dim_id)
            d = H.dimensionNumber(dim_id);
            head = H(d);
        end
        function [dim, dim_id] = non_singleton_dim(H)
            dim = find([H.n] > 1);
            dim_id = [H(dim).dim_id];
        end
        function dim_id = non_singleton_dim_id(H)
            dim_id = [H([H.n] > 1).dim_id];
        end
    end
    
end
