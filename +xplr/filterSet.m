classdef filterSet < hgsetget
    
    properties (SetAccess='private')
        linkkey
%     end
%     properties (Access='private')
        registry    % registered filter&point filters
        combo       % list combo displaying the registered filters
        zregistry   % registered zoomfilters
        zcregistry  % registered zoomcentral
    end
    
    % Constructor should be called only by xplr.bank.getFilterSet
    methods
        function FS = filterSet(key)
            FS.linkkey = key;
            FS.registry = xplr.bankRegistry;
            FS.zregistry = xplr.bankRegistry;
            FS.zcregistry = xplr.bankRegistry;
        end
    end
    
    % Filter registration
    methods
        function F = getFilter(FS,header, doshow, varargin)
            % function F = getFilter(FS,header, doshow [,user])
            hID = getID(header);
            F = FS.registry.getValue(hID,varargin{:});
            if isempty(F)
                F = xplr.filterAndPoint(header);
                % viewcontrol object C will be registered as a user of the filter
                FS.addFilter(F, doshow, varargin{:});
            end            
        end
        function F = getZoomFilter(FS,header,varargin)
            % function F = getZoomFilter(FS,header[,user])
            hID = getID(header);
            % search a corresponding filter
            F = FS.zregistry.getValue(hID,varargin{:});
            if ~isempty(F)
                return
            else
                % if the filter does not exist then create it and register it
                F = xplr.zoomfilter(header);
                FS.addZoomFilter(F,varargin{:});
                type = header.get.type();
                if strcmp(type, 'measure')
                   key = fn_hash({ header.label, header.get.unit }, 'num');
                   % search a corresponding zoomCentral in the filterSet
                   % return existing zoomCentral if exist
                   zoomCentral = FS.zcregistry.getValue(key,varargin{:});
                   if isempty(zoomCentral)
                       % create new zoomcentral
                       zoomCentral = xplr.zoomcentral(header.label,header.unit);
                       FS.zcregistry.register(key,zoomCentral,varargin{:});
                   end
                zoomCentral.connectZoomFilter(F);
                end
            end
        end
        function addFilter(FS,F,doshow,varargin)
            % function addFilter(FS,F[,user])
            hID = getID(F.headerin);
            F.linkkey = FS.linkkey;
            FS.registry.register(hID,F,varargin{:})
            if doshow 
                if F.ndin > 1
                    disp 'cannot display list for ND filter'
                else
                    FS.showList(F)
                end
            end
        end
        function addZoomFilter(FS,F,varargin)
            % function addZoomFilter(FS,F[,user])
            hID = getID(F.headerin);
            F.linkkey = FS.linkkey;
            FS.zregistry.register(hID,F,varargin{:})
        end
        function removeFilter(FS,F,varargin)
            % function removeFilter(FS,F[,user])
            hID = getID(F.headerin);
            removed = FS.registry.unregister(hID,varargin{:});
            if removed && ~isempty(FS.combo) && isvalid(FS.combo)
                FS.combo.removeList(F)
            end
        end
        function removeZoomFilter(FS,F,varargin)
            % function removeZoomFilter(FS,F[,user])
            if ~isvalid(FS)
                % can happen for example upon a 'clear all'
                return
            end
            hID = getID(F.headerin);
            FS.zregistry.unregister(hID,varargin{:});
        end
        function clear(FS)
            FS.registry.clear()
        end
    end
    
    % List display
    methods
        function showList(FS,F)
            % Create list combo?
            if isempty(FS.combo) || ~isvalid(FS.combo)
                FS.combo = xplr.listcombo([],FS.linkkey);
                % no need to delete the listener upon filterSet deletion: filterSet are supposed never to be deleted
                connectlistener(FS.combo,FS,'Empty',@(u,e)set(FS,'combo',[])); 
            end
            
            % Display filter if it is 1D
            if F.ndin == 1
                FS.combo.showList(F)
                figure(FS.combo.container.hobj)
            end
        end
    end
    
end