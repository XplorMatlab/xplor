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
    methods (Access='private')
        function datachangeSmart(S,dim,flag,ind)
            % Replace filter (as in slicer.replaceFilterDim)
            idx = find([S.filters.dim]==dim);
            oldfilt = S.filters(idx).obj;
            newfilt = autoZoomFilter(S,oldfilt.linkkey,dim); 
            % TODO: really create a new filter when flag is only 'chg'??
            % (to which slider needs to be reconnected, etc.)
            if strcmp(flag,'chg')
                % keep current zoom
                setZoom(newfilt,oldfilt.zoom,oldfilt.bin)
            end
            deleteValid(oldfilt,S.listeners.filters(idx))
            S.filters(idx).obj = newfilt;
            S.listeners.filters(idx) = addlistener(newfilt,'ChangedOperation',@(u,e)filterchange(S,dim,e));
            
            % Instead of reslicing everything, propagate changes in data in
            % a smart fashion
            switch flag
                case 'all'
                    % all data has changed, re-start slicing from zero
                    S.slicingchain(:) = [];
                    doslice(S,'data',flag,dim,ind)
                case {'remove' 'perm'}
                    % easy
                    for k = 1:length(S.slicingchain)
                        s = S.slicingchain(k);
                        s.res.updateData(flag,dim,ind) % last element will be the slice
                    end
                case {'new' 'chg' 'chg&new' 'chg&rm'}
                    % do a new slicing for a subset of the data

                    % get relevant subset of data
                    switch flag
                        case {'new' 'chg'}
                            ind1 = ind;
                        case 'chg&new'
                            ind1 = [ind{:}];
                        case 'chg&rm'
                            ind1 = ind{1};
                    end
                    headd = S.data.header(dim);
                    subs = substruct('()',repmat({':'},1,S.data.nd));
                    subs.subs{dim} = ind1;
                    datasub = subsref(S.data.data,subs);
                    
                    % slice it and update the slicing chain
                    prevres = S.data;
                    for k = 1:length(S.slicingchain)
                        s = S.slicingchain(k);
                        dimk = S.filters(k).dim;
                        filtk = S.filters(k).obj;
                        if ~strcmp(filtk.zoom,':') || filtk.bin~=1
                            % apply zoom in this dimension
                            if dim==dimk
                                if ~strcmp(flag,'chg'), error 'this is supposed to happen only with flag ''chg''!', end
                                zoomed = filtk.operation(prevres,dimk);
                                % indices of data change might not be the
                                % same anymore!!
                                ind = filtk.orig2zoomed(ind);
                                if ~any(ind), break, end % parts of data that have changed are not visible, we can stop here!
                                ind(~ind) = [];
                                headd = filtk.headerout;
                                s.res.updateData(flag,dim,ind,zoomed.data,headd) % last element is the slice
                            else
                                datasub = filtk.slicing(datasub,dimk);
                                s.res.updateData(flag,dim,ind,datasub,headd) % last element is the slice
                            end
                        else
                            s.res.updateData(flag,dim,ind,datasub,headd) % last element is the slice
                        end
                        prevres = s.res;
                    end
                otherwise
                    error('datachangeSmart cannot be called with flag ''%s''',flag)
            end
        end
    end
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
                        
                        % except for the case strcmp(e.flag,'chg'),
                        % dataChangeSmart does nothing more than
                        % slicer.doslice
                        %datachangeSmart(S,e.dim,e.flag,e.ind)

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