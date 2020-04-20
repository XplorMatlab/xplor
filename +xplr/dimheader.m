classdef dimheader < xplr.header
    % function H = dimheader(header arguments...)
    % function H = dimheader(header object[, dimID])
    %---
    % The DIMHEADER class inherits from HEADER and adds the 'dimID'
    % property, which is a random number serving to uniquely identify
    % different dimensions inside an xdata object, even if their headers
    % are rigorously the same.
    %
    % See also xplr.header
    
    properties
        dimID
    end
    
    methods
        function H = dimheader(varargin)
            % header part
            docopy = (nargin>=1 && isa(varargin{1},'xplr.header'));
            if docopy
                arg = {};
            else
                % construct header with xplr.header syntax
                arg = varargin;
            end
            H = H@xplr.header(arg{:});
            if docopy
                % copy header
                H1 = varargin{1};
                if ~isscalar(H1)
                    % pre-allocate
                    H(numel(H1)) = xplr.dimheader();
                    H = reshape(H,size(H1));
                end
                H.copyin(H1);
            end
            
            % generate a unique dimension ID
            if docopy && nargin>=2
                [H.dimID] = dealc(varargin{2});
            else
                H.changeDimID()
            end
        end
        function changeDimID(H)
            % generate a new, unique dimension ID; this method will be
            % called in particular when a header will be duplicated (and
            % potentially modified) for a new usage
            for i = 1:length(H)
                H(i).dimID = rand;
            end
        end
        function disp(H)
            disp@xplr.header(H)
            fprintf('\b    dimID: %.4f\n\n', H.dimID);
        end
        function [dim, dimID] = dimensionNumberAndID(H,d)
            % function [dim, dimID] = dimensionNumberAndID(H,d)
            %---
            % Convert any of dimension numbers, identifiers or labels
            % to both dimension numbers and identifiers.
            % Returns 0 for unfound dimensions.
            
            % Special cases
            if isempty(d)
                if iscell(d)
                    [dim, dimID] = deal(cell(1,0));
                else
                    [dim, dimID] = deal(zeros(1,0));
                end
                return
            elseif ischar(d) || (iscell(d) && ischar(d{1}))
                % First convert labels to dimension numbers
                if ~iscell(d), d = {d}; end
                n = length(d);
                dim = zeros(1,n);
                labels = {H.label};
                for i = 1:n
                    dim_i = fn_find(d{i},labels,'first');
                    if ~isempty(dim_i)
                        dim(i) = dim_i;
                    end
                end
                d = dim;
            elseif iscell(d)
                % Multiple outputs
                n = length(d);
                [dim, dimID] = deal(cell(1,n));
                for i = 1:n
                    [dim{i}, dimID{i}] = H.dimensionNumberAndID(d{i}); 
                end
                return
            end
            
            % Convert between dimension numbers and identifiers
            if d(1)<1
                % identifier -> number
                dimID = d;
                n = length(dimID);
                dim = zeros(1,n);
                for i =1:n
                    dim_i = find([H.dimID]==dimID(i),1,'first');
                    if ~isempty(dim_i)
                        dim(i) = dim_i;
                    end
                end
            else
                % number -> identifier
                dim = d;
                dimID = [H(dim).dimID];
                if length(dimID) < length(d)
                    [dim, dimID] = deal([]);
                    return
                end
            end
        end
        function dimID = dimensionID(H,d)
            [~, dimID] = H.dimensionNumberAndID(d);
        end
        function dim = dimensionNumber(H,d)
            [dim, ~] = H.dimensionNumberAndID(d);
        end
        function label = dimensionLabel(H,d)
            [dim, ~] = H.dimensionNumberAndID(d);
            if isempty(dim)
                label = [];
            else
                label = H(dim).label;
            end
        end
        function head = headerByID(H,dimID)
            d = H.dimensionNumber(dimID);
            head = H(d);
        end
        function [dim, dimID] = non_singleton_dim(H)
            dim = find([H.n]>1);
            dimID = [H(dim).dimID];
        end
        function dimID = non_singleton_dimID(H)
            dimID = [H([H.n]>1).dimID];
        end
    end
    
end