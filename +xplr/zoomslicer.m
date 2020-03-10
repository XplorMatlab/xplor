classdef zoomslicer < xplr.slicer
    % function S = zoomslicer(data)
    % ---
    % The zoomslicer class is a specialized version of the slicer class:
    % filters are not set by the user, but zoom filters are automatically
    % assigned to every dimension.
    
    properties (SetAccess='private')
        D           % parent 'display' object
        defaultlinkkey = 1;
    end
    
    events
        ChangedZoom % transit changes in the observed zoom filters
    end
    
    % Constructor
    methods 
        function S = zoomslicer(data,D)
            % Construct slicer object
            S = S@xplr.slicer(data);
            S.D = D;
            
            % Automatically create zoom filters for all dimensions
            Z = autoZoomFilter(S, S.defaultlinkkey);
            dim = 1:length(S.data.header);
            S.addFilter(dim,Z)
            
        end
        function Z = autoZoomFilter(S,linkkey,dim)
            % function Z = autoZoomFilter(S,linkkey[,dim])
            %---
            % Create, connect and return zoom filter for specified (or all)
            % dimension(s) and for specified link key (get from / add to
            % the bank if appropriate).
            % But does not add / replace it in the zoomslicer: this needs
            % to be done by the calling function.
            
            % to do: check the bank for available filter instead of
            % creating systematically a new one
            if nargin<3, dim = 1:length(S.data.header); end
            
            for i=1:length(dim)
                d = dim(i);
                head = S.data.header(d);
                if linkkey ~= 0
                    Zi = xplr.bank.getZoomFilter(linkkey,head,S);
                else
                    Zi = xplr.zoomfilter(S.data.header(dim));
                end
                Z(i) = Zi;
                S.addListener(Z(i),'ChangedZoom',@(u,e)zoomchange(S,d,e));
            end
        end       
    end
    
    % Action upon data change differ from the parent slicer class
    methods
        function datachange(S,e)
            switch e.flag
                case 'global'
                    % remove all existing filters and create new ones
                    S.rmFilter(1:length(S.filters))
                    Z = autoZoomFilter(S,S.defaultlinkkey);
                    dim = 1:length(S.data.header);
                    S.addFilter(dim,Z)
                case 'chgdim'
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    newfilt = autoZoomFilter(S,S.defaultlinkkey,e.dim);
                    for i=1:length(e.dim)
                        replaceFilterDim(S,e.dim(i),newfilt(i),false)
                    end
                    doslice(S,'data','chgdim',e.dim)
                case {'all' 'new' 'remove' 'chg&new' 'chg&rm' 'perm' 'chg'}
                    % In this case we do not use methods of the parent
                    % 'slicer' class, but do all the update here in a
                    % smarter way...
                    dim = e.dim;
                    idx = ([S.filters.dim]==dim);
                    curfilt = S.filters(idx).obj;
                    if strcmp(e.flag,'chg') || (strcmp(curfilt.zoom,':') && curfilt.bin==1)
                        % 'smart' filter replacement: partial update will
                        % be possible
                        
                        % replace filter (as in slicer.replaceFilterDim)
                        S.disconnect(curfilt)
                        newfilt = autoZoomFilter(S,curfilt.linkkey,dim);
                        % TODO: really create a new filter when flag is only 'chg'??
                        % (to which slider needs to be reconnected, etc.)
                        if strcmp(e.flag,'chg')
                            % keep current zoom
                            setZoom(newfilt,curfilt.zoom,curfilt.bin)
                        end
                        S.filters(idx).obj = newfilt;
                        S.addListener(newfilt,'ChangedOperation',@(u,e)filterchange(S,dim,e));

                        % smart update
                        if strcmp(e.flag,'all'), ind = []; else ind = e.ind; end
                        doslice(S,'data',e.flag,e.dim,ind)
                    else
                        % full update
                        S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                        newfilt = autoZoomFilter(S,curfilt.linkkey,e.dim);
                        replaceFilterDim(S,e.dim,newfilt)
                    end
                case 'chgdata'
                    % no change in the input header
                    S.slicingchain(:) = [];
                    doslice(S,'data','chgdata')
                otherwise
                    error 'not implemented yet'
            end
        end
    end
    
    % Transit zoom changes
    methods
        function setZoom(S,dim,newzoom)
            c = S.disableConnection(S.filters(dim));
            chgnout = false;
            for i=1:length(dim)
                Z = S.filters(dim(i)).obj;
                curnout = Z.szout;
                Z.setZoom(newzoom(:,i))
                chgnout = chgnout || (Z.szout~=curnout);
            end
            notify(S,'ChangedZoom',xplr.eventinfo('zoom',chgnout,dim)) % to do: check whether chgnout is true or false...
        end
        function zoomchange(S,dim,e)
            notify(S,'ChangedZoom',xplr.eventinfo('zoom',e.chgnout,dim))
        end
    end
    
 

    
    % Automatic unregistration from the bank upon disconnection
    methods
        function disconnect(S,F)
            disconnect@xplr.graphnode(S,F)
            if isa(F,'xplr.zoomfilter') && F.linkkey ~= 0
                xplr.bank.unregisterZoomFilter(F,S)
            end
        end
    end
    
    % Replace a filter by a new one with another linkkey
    methods
        function changeKey(S,dim,key)
            S.replaceFilterDim(dim,S.autoZoomFilter(key,dim),1);
        end
    end
end