classdef ZoomSlicer < xplr.Slicer
    % function S = zoomslicer(data)
    % ---
    % The zoomslicer class is a specialized version of the slicer class:
    % filters are not set by the user, but zoom filters are automatically
    % assigned to every dimension. So there is always as many filters as
    % dimension (even singleton ones).

    % A difficulty in the code of the zoomslicer is that dimensions are
    % sometimes identified by their number (dim), or by their identifier
    % (dim_id): developpers should not get confused between the two of them!
    % Some methods accept both identification methods (dimension numbers
    % can easily be distinguished from dimension identifiers because the
    % former are >=1, the latter are <1).
    
    properties (SetAccess='private')
        D           % parent 'display' object
        default_link_key = 1;
    end
    
    events
        ChangedZoom % transit changes in the observed zoom filters
    end
    
    % Constructor
    methods 
        function S = ZoomSlicer(data, D)
            % Construct slicer object
            S = S@xplr.Slicer([], data);
            S.D = D;
            
            % Automatically create zoom filters for all dimensions
            Z = auto_zoom_filter(S, S.default_link_key);
            dim = 1:length(S.data.header);
            S.add_filter(dim, Z)
            
        end
    end
    
    % Automatic creation of zoom filters (contrary to the Slicer parent
    % class that receives only filters defined externally, the ZoomSlicer
    % object is the one that creates and manages its filters)
    methods
        function Z = auto_zoom_filter(S, link_key, dim)
            % function Z = auto_zoom_filter(S,link_key[,dim])
            %---
            % Create, connect and return zoom filter for specified (or all)
            % dimension(s) and for specified link key (get from / add to
            % the bank if appropriate).
            % But does not add / replace it in the zoomslicer: this needs
            % to be done by the calling function.
            
            % to do: check the bank for available filter instead of
            % creating systematically a new one
            if nargin < 3, dim = 1:length(S.data.header); end
            
            for i=1:length(dim)
                d = dim(i);
                head = S.data.header(d);
                if link_key ~= 0
                    zi = xplr.Bank.get_zoom_filter(link_key, head, S);
                else
                    zi = xplr.ZoomFilter(S.data.header(dim));
                end
                Z(i) = zi;
                S.add_listener(Z(i), 'ChangedOperation', @(u,e)zoom_filter_change(S,d,e));
            end
        end       
    end
    
    % Action upon data change differ from the parent slicer class
    methods
        function data_change(S,e)
            dim_id = e.dim;
            dim = S.data.dimension_number(dim_id);
            if length(S.filters) ~= S.data.nd || ~all(isvalid([S.filters.obj]))
                % something is wrong, use 'global' flag
                e.flag = 'global';
            end
            switch e.flag
                case 'global'
                    % remove all existing filters and create new ones
                    S.slicing_chain(:) = []; % data has changed, all previous slicing steps became invalid
                    S.rm_filter(1:length(S.filters), false) % no need to reslice at this stage, there will be a reslice below
                    Z = auto_zoom_filter(S, S.default_link_key);
                    dim = 1:length(S.data.header);
                    S.add_filter(dim, Z)
                case {'chg_dim', 'all'}
                    % replace filters in modified dimensions
                    S.slicing_chain(:) = []; % data and dimensions have changed, all previous slicing steps became invalid
                    Z = auto_zoom_filter(S, S.default_link_key, dim);
                    S.replace_filter(dim, dim, Z)
                case {'new', 'remove', 'chg&new', 'chg&rm', 'perm', 'chg'}
                    % In this case we do not use methods of the parent
                    % 'slicer' class, but do all the update here in a
                    % smarter way...
                    cur_filt = S.filters(dim).obj;
                    if strcmp(cur_filt.zoom, ':') && cur_filt.bin == 1
                        % the current filter has no effect (no zooming, no
                        % binning), so we can propagate the smart update
                        % information from the slice to the zslice, where
                        % the concerned dimension will share the same dim_id
                        % (xplr.dataOperand.changedimension_id)
                        
                        % replace filter (as in slicer.replace_filter_dim)
                        cur_key = cur_filt.link_key;
                        S.disconnect(cur_filt) % cur_filt will be deleted if it is not used elsewhere
                        new_filt = auto_zoom_filter(S, cur_key, dim);
                        S.filters(dim).obj = new_filt;
                        S.add_listener(new_filt, 'ChangedOperation', @(u,e)filter_change(S,new_filt,dim,e));

                        % smart update: note that this will call
                        % maybe we need: S.slicing_chain(dim:end) = [];
                        if strcmp(e.flag, 'all'), ind = []; else, ind = e.ind; end
                        do_slice(S, 'data', e.flag,dim,ind)
                    else
                        % full update
                        S.slicing_chain(:) = []; % data has changed, all previous slicing steps became invalid
                        new_filt = auto_zoom_filter(S, cur_filt.link_key, dim);
                        replace_filter_dim(S, dim, new_filt)
                    end
                case 'chg_data'
                    % no change in the input header
                    do_slice(S, 'data', 'chg_data')
                otherwise
                    error 'not implemented yet'
            end
        end
    end
    
    % Transit zoom changes
    methods
        function set_zoom(S, dim, new_zoom)
            c = S.disable_connection(S.filters(dim));
            chg_n_out = false;
            for i=1:length(dim)
                Z = S.filters(dim(i)).obj;
                cur_n_out = Z.sz_out;
                Z.set_zoom(new_zoom(:, i))
                chg_n_out = chg_n_out || (Z.sz_out ~= cur_n_out);
            end
            notify(S, 'ChangedZoom', xplr.EventInfo('zoom', chg_n_out,dim)) % to do: check whether chg_n_out is true or false...
        end
        function zoom_filter_change(S, dim, e)
            if strcmp(e.type, 'zoom')
                notify(S, 'ChangedZoom', xplr.EventInfo('zoom', e.chg_n_out, dim))
            end
        end
    end
    
 

    
    % Automatic unregistration from the bank upon disconnection
    methods
        function disconnect(S, F)
            % Multiple filters
            if ~isscalar(F)
                for i = 1:length(F)
                    disconnect(S, F(i))
                end
                return
            end
            
            % Disconnect one filter
            disconnect@xplr.GraphNode(S, F)
            if isa(F, 'xplr.ZoomFilter') && isvalid(F) && F.link_key ~= 0
                xplr.Bank.unregister_filter(F, S)
            end
        end
    end
    
    % Replace a filter by a new one with another link_key
    methods
        function change_key(S, dim, key)
            S.replace_filter_dim(dim, S.auto_zoom_filter(key,dim));
        end
    end
end
