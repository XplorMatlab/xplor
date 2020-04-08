classdef zoomslicer < xplr.slicer
    % function S = zoomslicer(data)
    % ---
    % The zoomslicer class is a specialized version of the slicer class:
    % filters are not set by the user, but zoom filters are automatically
    % assigned to every dimension. So there is always as many filters as
    % dimension (even singleton ones).

    % A difficulty in the code of the zoomslicer is that dimensions are
    % sometimes identified by their number (dim), or by their identifier
    % (dimID): developpers should not get confused between the two of them!
    % Some methods accept both identification methods (dimension numbers
    % can easily be distinguished from dimension identifiers because the
    % former are >=1, the latter are <1).
    
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
            S = S@xplr.slicer([],data);
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
                S.addListener(Z(i),'ChangedOperation',@(u,e)zoomfilterchange(S,d,e));
            end
        end       
    end
    
    % Action upon data change differ from the parent slicer class
    methods
        function datachange(S,e)
            dimID = e.dim;
            dim = S.data.dimensionNumber(dimID);
            if length(S.filters) ~= S.data.nd || ~all(isvalid([S.filters.obj]))
                % something is wrong, use 'global' flag
                e.flag = 'global';
            end
            switch e.flag
                case 'global'
                    % remove all existing filters and create new ones
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    S.rmFilter(1:length(S.filters),false) % no need to reslice at this stage, there will be a reslice below
                    Z = autoZoomFilter(S,S.defaultlinkkey);
                    dim = 1:length(S.data.header);
                    S.addFilter(dim,Z)
                case {'chgdim' 'all'}
                    % replace filters in modified dimensions
                    S.slicingchain(:) = []; % data and dimensions have changed, all previous slicing steps became invalid
                    Z = autoZoomFilter(S,S.defaultlinkkey,dim);
                    S.replaceFilter(dim,dim,Z)
                case {'new' 'remove' 'chg&new' 'chg&rm' 'perm' 'chg'}
                    % In this case we do not use methods of the parent
                    % 'slicer' class, but do all the update here in a
                    % smarter way...
                    curfilt = S.filters(dim).obj;
                    if strcmp(curfilt.zoom,':') && curfilt.bin==1
                        % the current filter has no effect (no zooming, no
                        % binning), so we can propagate the smart update
                        % information from the slice to the zslice, where
                        % the concerned dimension will share the same dimID
                        % (xplr.dataOperand.changedimensionID)
                        
                        % replace filter (as in slicer.replaceFilterDim)
                        curkey = curfilt.linkkey;
                        S.disconnect(curfilt) % curfilt will be deleted if it is not used elsewhere
                        newfilt = autoZoomFilter(S,curkey,dim);
                        S.filters(dim).obj = newfilt;
                        S.addListener(newfilt,'ChangedOperation',@(u,e)filterchange(S,newfilt,dim,e));

                        % smart update: note that this will call
                        % maybe we need: S.slicingchain(dim:end) = [];
                        if strcmp(e.flag,'all'), ind = []; else, ind = e.ind; end
                        doslice(S,'data',e.flag,dim,ind)
                    else
                        % full update
                        S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                        newfilt = autoZoomFilter(S,curfilt.linkkey,dim);
                        replaceFilterDim(S,dim,newfilt)
                    end
                case 'chgdata'
                    % no change in the input header
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
        function zoomfilterchange(S,dim,e)
            if strcmp(e.type,'zoom')
                notify(S,'ChangedZoom',xplr.eventinfo('zoom',e.chgnout,dim))
            end
        end
    end
    
 

    
    % Automatic unregistration from the bank upon disconnection
    methods
        function disconnect(S,F)
            % Multiple filters
            if ~isscalar(F)
                for i = 1:length(F)
                    disconnect(S,F(i))
                end
                return
            end
            
            % Disconnect one filter
            disconnect@xplr.graphnode(S,F)
            if isa(F,'xplr.zoomfilter') && isvalid(F) && F.linkkey ~= 0
                xplr.bank.unregisterFilter(F,S)
            end
        end
    end
    
    % Replace a filter by a new one with another linkkey
    methods
        function changeKey(S,dim,key)
            S.replaceFilterDim(dim,S.autoZoomFilter(key,dim));
        end
    end
end