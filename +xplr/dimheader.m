classdef dimheader < xplr.header
    % function H = dimheader(header arguments...)
    % function H = dimheader(header object)
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
            docopy = (nargin==1 && isa(varargin{1},'xplr.header'));
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
            H.changeDimID()
        end
        function changeDimID(H)
            % generate a new, unique dimension ID; this method will be
            % called in particular when a header will be duplicated (and
            % potentially modified) for a new usage
            for i = 1:length(H)
                H(i).dimID = rand;
            end
        end
    end
    
end