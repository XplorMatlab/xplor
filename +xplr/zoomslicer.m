classdef zoomslicer < xplr.slicer
    % function S = zoomslicer(data)
    %---
    % The zoomslicer class is a specialized version of the slicer class:
    % filters are not set by the user, but zoom filters are automatically
    % assigned to every dimension.
    
    properties (SetAccess='private')
        D           % parent 'display' object
        defaullinkkey = 1;
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
            autoZoomFilter(S, S.defaultlinkkey)
        end
        function Z = autoZoomFilter(S,linkkey,dim)
            % function Z = autoZoomFilter(S,linkkey[,dim])
            %---
            % Create zoom filter for specified (or all) dimension(s). Add
            % the filter(s) only if there is no output requested.
            
            % to do: check the bank for available filter instead of
            % creating systematically a new one
            if nargin<2, dim = 1:length(S.data.header); end
            Z = xplr.zoomfilter(S,S.data.header(dim));
            
            if ~isfield(S.listeners,'zoom'), S.listeners.zoom = event.listener.empty(1,0); end
            for i=1:length(dim)
                d = dim(i);
                if length(S.listeners.zoom)>=d, deleteValid(S.listeners.zoom(d)), end
                S.listeners.zoom(d) = addlistener(Z(i),'ChangedZoom',@(u,e)zoomchange(S,d,e));
            end
            if nargout==0
                S.addFilter(dim,Z);
                % register the zoom filter in bank.filterSet(key).zregistry
                xplr.bank.addZoomFilter(1, Z, S.D.V.C)
                clear Z
            end
            
            
        end
        
    end
    
    % Action upon data change differ from the parent slicer class
    methods
        function datachange(S,e)
            switch e.flag
                case 'global'
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    S.filters(:) = [];
                    autoZoomFilter(S,S.defaultlinkkey)
                case 'chgdim'
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    newfilt = autoZoomFilter(S,S.defaultlinkkey,e.dim);
                    for i=1:length(e.dim)
                        replaceFilterDim(S,e.dim(i),newfilt(i),false)
                    end
                    doslice(S,'data','chgdim',e.dim)
                case {'all' 'new' 'remove' 'chg&new' 'chg&rm' 'perm' 'chg'}
                    dim = e.dim;
                    idx = ([S.filters.dim]==dim);
                    curfilt = S.filters(idx).obj;
                    if strcmp(e.flag,'chg') || (strcmp(curfilt.zoom,':') && curfilt.bin==1)
                        % 'smart' filter replacement: partial update will
                        % be possible
                        
                        % replace filter (as in slicer.replaceFilterDim)
                        newfilt = autoZoomFilter(S,curfilt.linkkey,dim);
                        % TODO: really create a new filter when flag is only 'chg'??
                        % (to which slider needs to be reconnected, etc.)
                        if strcmp(e.flag,'chg')
                            % keep current zoom
                            setZoom(newfilt,curfilt.zoom,curfilt.bin)
                        end
                        deleteValid(curfilt,S.listeners.filters(idx))
                        S.filters(idx).obj = newfilt;
                        S.listeners.filters(idx) = addlistener(newfilt,'ChangedOperation',@(u,e)filterchange(S,dim,e));

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
            c = disableListener(S.listeners.zoom(dim)); %#ok<NASGU>
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
    
end