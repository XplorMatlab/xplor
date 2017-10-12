classdef filterSet < hgsetget
    
    properties (SetAccess='private')
        linkkey
%     end
%     properties (Access='private')
        registry
        combo
    end
    
    % Constructor should be called only by xplr.bank.getFilterSet
    methods
        function S = filterSet(key)
            S.linkkey = key;
            S.registry = xplr.bankRegistry;
        end
    end
    
    % Filter registration
    methods
        function F = getFilter(S,header,varargin)
            % function F = getFilter(S,header[,user])
            hID = getID(header);
            F = S.registry.getValue(hID,varargin{:});
        end
        function addFilter(S,F,varargin)
            % function addFilter(S,F[,user])
            hID = getID(F.headerin);
            F.linkkey = S.linkkey;
            S.registry.register(hID,F,varargin{:})
            % Automatic unregister upon user's deletion
            if nargin>=3
                user = varargin{1};
                fn_deletefcn(user,@(u,e)removeFilter(S,F,user))
            end
        end
        function removeFilter(S,F,varargin)
            % function removeFilter(S,F[,user])
            hID = getID(F.headerin);
            removed = S.registry.unregister(hID,varargin{:});
            if removed
                S.combo.removeList(F)
            end
        end
    end
    
    % List display
    methods
        function showList(S,F)
            % Create list combo?
            if isempty(S.combo)
                S.combo = xplr.listcombo([],S.linkkey);
                connectlistener(S.combo,S,'Empty',@(u,e)set(S,'combo',[])); % no need to delete the listener upon filterSet deletion: filterSet are supposed never to be deleted
            end
            
            % Display filter
            S.combo.showList(F)
        end
    end
    
end