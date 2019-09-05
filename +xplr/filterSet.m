classdef filterSet < hgsetget
    
    properties (SetAccess='private')
        linkkey
%     end
%     properties (Access='private')
        registry    % registered filters
        combo       % list combo displaying the registered filters
        zregistry   % registered zoomfilters
    end
    
    % Constructor should be called only by xplr.bank.getFilterSet
    methods
        function FS = filterSet(key)
            FS.linkkey = key;
            FS.registry = xplr.bankRegistry;
            FS.zregistry = xplr.bankRegistry;
        end
    end
    
    % Filter registration
    methods
        function F = getFilter(FS,header,varargin)
            % function F = getFilter(FS,header[,user])
            hID = getID(header);
            F = FS.registry.getValue(hID,varargin{:});
            % Automatic unregister upon user's deletion
            if ~isempty(F) && nargin>=3
                user = varargin{1};
                fn_deletefcn(user,@(u,e)removeFilter(FS,F,user))
            end
        end
        function F = getZoomFilter(FS,header,varargin)
            % function F = getZoomFilter(FS,header[,user])
            hID = getID(header);
            F = FS.zregistry.getValue(hID,varargin{:});
            % Automatic unregister upon user's deletion
            if ~isempty(F) && nargin>=3
                user = varargin{1};
                fn_deletefcn(user,@(u,e)removeZoomFilter(FS,F,user))
            end
        end
        function addFilter(FS,F,varargin)
            % function addFilter(FS,F[,user])
            hID = getID(F.headerin);
            F.linkkey = FS.linkkey;
            FS.registry.register(hID,F,varargin{:})
            FS.showList(F)
            % Automatic unregister upon user's deletion
            if nargin>=3
                user = varargin{1};
                fn_deletefcn(user,@(u,e)removeFilter(FS,F,user))
            end
        end
        function addZoomFilter(FS,F,varargin)
            % function addZoomFilter(FS,F[,user])
            hID = getID(F.headerin);
            F.linkkey = FS.linkkey;
            FS.zregistry.register(hID,F,varargin{:})
            % Automatic unregister upon user's deletion
            if nargin>=3
                user = varargin{1};
                fn_deletefcn(user,@(u,e)removeZoomFilter(FS,F,user))
            end
        end
        function removeFilter(FS,F,varargin)
            % function removeFilter(FS,F[,user])
            hID = getID(F.headerin);
            removed = FS.registry.unregister(hID,varargin{:});
            if removed && ~isempty(FS.combo)
                FS.combo.removeList(F)
            end
        end
        function removeZoomFilter(FS,F,varargin)
            % function removeZoomFilter(FS,F[,user])
            hID = getID(F.headerin);
            removed = FS.zregistry.unregister(hID,varargin{:});
            if removed && ~isempty(FS.combo)
                FS.combo.removeList(F)
            end
        end
        function clear(FS)
            FS.registry.clear()
        end
    end
    
    % List display
    methods
        function showList(FS,F)
            % Create list combo?
            if isempty(FS.combo)
                FS.combo = xplr.listcombo([],FS.linkkey);
                connectlistener(FS.combo,FS,'Empty',@(u,e)set(FS,'combo',[])); % no need to delete the listener upon filterSet deletion: filterSet are supposed never to be deleted
            end
            
            % Display filter
            FS.combo.showList(F)
            figure(FS.combo.container.hobj)
        end
    end
    
end