classdef Slicer < xplr.GraphNode
    % function S = slicer(data[,dims,filters])
    % Compute and automatically update operations applied to an initial
    % 'data', resulting in an output data called 'slice'.
    % In addition to filters, points define positions where to extract the
    % data in some dimensions in the event where there is no filter in this
    % dimension, or the filter is empty.
    
    % In general, data dimensions are identified by their unique identifier
    % (dim_id) rather than by their number, but some methods accept both
    % identification methods (dimension numbers can easily be distinguished
    % from dimension identifiers because the former are >=1, the latter are
    % <1 ; conversion of either number, identifier or label to identifier
    % is performed by the xplr.xdata.dimension_id method).
    
    properties 
        V
        data
        slice
        filters = struct('active', [], 'dim_id', cell(1,0), 'obj', []);
    end
    properties (SetAccess='protected')
        slicing_chain = xplr.XData.empty(1, 0); % intermediary and last steps of slicing, length is same as S.active_filters
        pend_ing_rm_filter = false; % remember when some dimensions were removed but reslice did not occur yet
    end
    properties (Dependent, SetAccess='private')
        nd_data
        nd_slice
        n_active_filt
        active_filters
    end
    
    % Constructor, destructor, basic access and get/set dependent
    methods
        function S = Slicer(V, data, dim_id, filters)
            % link to parent view
            S.V = V;
            xplr.debug_info('TODO', 'can a slicer exist without being aware of its possessing view?')
            % set data
            S.data = data;
            S.add_listener(data, 'changed_data', @(u,e)data_change(S,e));
            % without any filter, slice is identical data
            S.slice = data.copy();

            % set filters
            if nargin >= 3 && ~isempty(filters)
                add_filter(S, dim_id, filters)
            end
        end
        function n = get.nd_data(S)
            n = S.data.nd;
        end
        function n = get.nd_slice(S)
            n = S.slice.nd;
        end
        function n = get.n_active_filt(S)
            n = sum([S.filters.active]);
        end
        function filters = get.active_filters(S)
            filters = S.filters([S.filters.active]);
        end
    end

    % Changing filters
    methods
        function insert_filter(S, id_x, dim_id, new_filt, active)
            % function insert_filter(S,id_x,dim_id,new_filt[active])
            dim_id = S.data.dimension_id(dim_id);
            if nargin < 5, active = true; end
            
            % check filter
            if ~isa(new_filt, 'xplr.DataOperand'), error 'filter must be a dataOperand object', end
            
            % check dimensions and headers
            if ~iscell(dim_id)
                if isscalar(new_filt), dim_id = {dim_id}; else dim_id = num2cell(dim_id); end
            end
            nadd = length(dim_id);
            dim_id_add = [dim_id{:}];
            if length(unique(dim_id_add)) < length(dim_id_add)
                error 'same dimension appears in multiple filters'
            elseif any(ismember(dim_id_add, [S.filters.dim_id]))
                error 'some filter is already applied in at least one of the specified dimensions'
            elseif ~isequal([new_filt.header_in], S.data.header_by_id(dim_id_add))
                error 'some filter header does not match data'
            end
            % check header
            % add filters
            if isempty(id_x), id_x = length(S.filters) + 1; end
            S.filters = [S.filters(1:id_x-1), struct('active',num2cell(active),'dim_id',dim_id,'obj',num2cell(new_filt)), S.filters(id_x:end)];
            % remove invalid part of the slicing chain
            nok = sum([S.filters(1:id_x-1).active]);
            S.slicing_chain(nok+1:end) = [];
            % install listeners
            for i=1:nadd
                S.add_listener(new_filt(i), 'changed_operation', @(u,e)filter_change(S, new_filt(i), dim_id{i}, e));
            end
            % update slice
            if ~S.pend_ing_rm_filter && isscalar(new_filt) && new_filt.nd_out == new_filt.nd_in
                do_slice(S,'slicer', 'chg_dim',dim_id_add)
            else
                do_slice(S, 'slicer', 'global')
                S.pend_ing_rm_filter = false;
            end
        end
        function add_filter(S, dim_id, new_filt, active)
            % function add_filter(S,dim_id,new_filt[,active])
            if nargin < 4, active = true; end
            insert_filter(S, [], dim_id, new_filt, active)
        end
        function rm_filter(S, id_x, do_slicing)
            % function rm_filter(S,id_x[,do_slicing])
            if nargin < 3, do_slicing = true; end
            filt_rm = [S.filters(id_x).obj];
            S.disconnect(filt_rm)
            % among the dimensions for which filters will be removed, which
            % are those where the filters were really active
            active = [S.filters(id_x).active];
            id_xactive = id_x(active);
            chg_dim_id = [S.filters(id_xactive).dim_id];
            % remove invalid part of slicing chain
            if ~isempty(id_xactive)
                nok = sum([S.filters(1:id_xactive(1)-1).active]);
                S.slicing_chain(nok+1:end) = [];
            end
            % remove the filters
            S.filters(id_x) = [];
            % update slice
            if isempty(id_xactive), return, end
            if do_slicing
                if all(isvalid(filt_rm)) && all([filt_rm(active).nd_out] == [filt_rm(active).nd_in])
                    do_slice(S, 'slicer','chg_dim',chg_dim_id)
                else
                    do_slice(S, 'slicer', 'global')
                end
            else
                % remember that reslice did not occur after removing filter
                S.pend_ing_rm_filter = true;
            end
        end
        function rm_filter_dim(S, dim_id, do_slicing)
            % function rm_filter_dim(S,dim_id[,do_slicing])
            dim_id = S.data.dimension_id(dim_id);
            if nargin < 3, do_slicing = true; end
            % remove filters
            id_xrm = S.get_filter_index(dim_id);
            if ~any(id_xrm)
                fprintf('there is already no filter in dimension %i! aborting rm_filter_dim\n', dim_id)
                return
            end
            rm_filter(S, id_xrm, do_slicing)
        end
        function replace_filter(S, id_xs, dim_ids, new_filt)
            % function replace_filter(S,id_xs,dim_ids,new_filt)
            
            % there can be several filters: dim_ids must be a cell array
            dim_ids = S.data.dimension_id(dim_ids);
            if ~iscell(dim_ids)
                if ~isscalar(new_filt), dim_ids = num2cell(dim_ids); else, dim_ids = {dim_ids}; end
            end    
            if any(diff([length(id_xs), length(dim_ids), length(new_filt)])), error 'length mismatch', end

            % replace filter(s)
            nd_out_changed = false;
            chg_dim = zeros(1, 0);
            for i = 1:length(dim_ids)
                dim_id = dim_ids{i};
                F = new_filt(i);
                id_x = id_xs(i);
                % replace filter
                prev_nd_out = S.filters(id_x).obj.nd_out;
                S.disconnect(S.filters(id_x).obj)
                S.filters(id_x).obj = F;
                S.filters(id_x).dim_id = dim_id;
                S.add_listener(F, 'changed_operation', @(u,e)filter_change(S, F, dim_id, e));
                % no need for update if the filter is not active
                if S.filters(id_x).active
                    if F.nd_out ~= prev_nd_out, nd_out_changed = true; end
                    chg_dim = [chg_dim, dim_id]; %#ok<AGROW>
                else
                    S.activate_connection(F, false)
                end
                % remove invalid part of slicing chain
                nok = sum([S.filters(1:id_x(1)-1).active]);
                S.slicing_chain(nok+1:end) = [];
            end
            
            % update slice
            if nd_out_changed
                do_slice(S, 'slicer', 'global')
            else
                do_slice(S, 'slicer', 'chg_dim', chg_dim)
            end
        end
        function replace_filter_dim(S, dim_ids, new_filt)
            % function replace_filter_dim(S,dims,new_filt)
            
            % there can be several filters: dim_ids must be a cell array
            dim_ids = S.data.dimension_id(dim_ids);
            if ~iscell(dim_ids)
                if ~isscalar(new_filt), dim_ids = num2cell(dim_ids); else, dim_ids = {dim_ids}; end
            end               

            % replace filter(s)
            nd_out_changed = false;
            chg_dim = zeros(1, 0);
            for i = 1:length(dim_ids)
                dim_id = dim_ids{i};
                F = new_filt(i);
                % replace filter
                id_x = fn_find(dim_id, {S.filters.dim_id});
                if isempty(id_x), error 'there is no filter in the specified dimension', end
                prev_nd_out = S.filters(id_x).obj.nd_out;
                S.disconnect(S.filters(id_x).obj)
                S.filters(id_x).obj = F;
                S.add_listener(F, 'changed_operation', @(u,e)filter_change(S, F, dim_id, e));
                % no need for update if the filter is not active
                if S.filters(id_x).active
                    if F.nd_out~=prev_nd_out, nd_out_changed = true; end
                    chg_dim = [chg_dim, dim_id]; %#ok<AGROW>
                else
                    S.activate_connection(F, false)
                end
                % remove invalid part of slicing chain
                nok = sum([S.filters(1:id_x(1)-1).active]);
                S.slicing_chain(nok+1:end) = [];
            end
            
            % update slice
            if nd_out_changed
                do_slice(S, 'slicer', 'global')
            else
                do_slice(S, 'slicer', 'chg_dim', chg_dim)
            end
        end
        function perm_filters(S, perm)
            % function perm_filters(S,perm)
            % modify the order in which filters are applied
            
            % check
            if ~isequal(sort(perm),1:length(S.filters))
                error 'argument is not a valid permutation'
            end
            id_xfirst = find(perm~=1:length(S.filters), 1, 'first');
            if isempty(id_xfirst)
                % no permutation at all!
                return
            end
            
            % permutation
            S.filters = S.filters(perm);
            
            % update slice: note that all output header will remain the
            % same, so best is to signal the change as a change in the data
            nok = sum([S.filters(1:id_xfirst-1).active]);
            S.slicing_chain(nok+1:end) = [];
            do_slice(S, 'slicer', 'chgdata')
        end
        function chg_filter_active(S, id_x, val)
            % function chg_filter_active(S,id_x,val)
            val = boolean(val);
            for i = id_x
                S.filters(i).active = val;
                S.activate_connection(S.filters(i).obj, val)
            end
            nok = sum([S.filters(1:min(id_x)-1).active]);
            S.slicing_chain(nok+1:end) = [];
            do_slice(S, 'slicer', 'global')
        end
        function applyPend_ing(S)
            if S.pend_ing_rm_filter
                S.do_slice('slicer', 'global')
                S.pend_ing_rm_filter = false;
            end
        end
    end

    % Get filter
    methods
        function F = get_filter_by_dim(S, dim_id)
            dim_id = S.data.dimension_id(dim_id);
            id_x = fn_find(dim_id, {S.filters.dim_id});
            F = [S.filters(id_x).obj];
        end
        function id_x = get_filter_index(S, F_or_dim)
            % function id_x = get_filter_index(S,[F|dim|dim_id])
            if isnumeric(F_or_dim)
                dim_id = S.data.dimension_id(F_or_dim);
                id_x = false(1, length(S.filters));
                for i=1:length(S.filters)
                    id_x(i) = any(ismember(S.filters(i).dim_id, dim_id));
                end
                id_x = find(id_x);
            else
                F = F_or_dim;
                if ~isscalar(F), error 'input filter must be scalar', end
                id_x = fn_find(F, [S.filters.obj], 'first');
            end
        end
    end
    
    % Slicing
    methods (Access='protected')
        function do_slice(S, chg_origin, chg_flag, chg_dim, chg_ind, chg_filter)
            % function do_slice(S,chg_origin,flag[,chg_dim[,ind[,filter]]])
            %---
            % Apply filters successively to get a slice out of the data.
            % "Smart updates" apply when the flag (for example 'new')
            % indicates that only a part of the data will need to be
            % re-calculated (for example with 'new', all the existing slice
            % will be preserved, but some additional part will be added
            % correspond_ing to the original increase in data or filter).

            if nargin < 4, chg_dim = []; end
            if nargin < 5, chg_ind = []; end
            if nargin < 6, chg_filter = []; end
            
            n_active_filters = sum([S.filters.active]);
            
            % 1) Smart update of the existing part of the slicing chain
            % dimension number
            [chg_dim, chg_dim_id] = S.data.dimension_number_and_id(chg_dim);
            % indices of the change
            switch chg_flag
                case {'new', 'chg'}
                    ind1 = chg_ind;
                case 'chg&new'
                    ind1 = [chg_ind{:}];
                case 'chg&rm'
                    ind1 = chg_ind{1};
                otherwise
                    ind1 = [];
            end
            % where to start, and get the relevant sub-part of the data
            % needed for smart update calculations
            switch chg_origin
                case 'data'
                    k_start = 1;
                    if fn_ismemberstr(chg_flag,{'all', 'chgdata', 'global', 'chg_dim'})
                        % smart update is not possible
                        S.slicing_chain(:) = [];
                    else
                        if ~isscalar(chg_dim), error 'programming', end
                        % original data has changed in dimension dim_id
                        head_dim = S.data.header(chg_dim);
                        subs = substruct('()', repmat({':'},1,S.data.nd));
                        subs.subs{chg_dim} = ind1;
                        % part of the data that will be used for new
                        % calculations, while preserving some of the
                        % existing calculations
                        data_sub = subsref(S.data.data, subs);
                        % the 'smart' update however will not be possible
                        % once there will be another filter applied in the
                        % same dimension dim
                        for k=1:length(S.slicing_chain)
                            filter_k = S.active_filters(k);
                            if filter_k.active && any(intersect(filter_k.dim_id, chg_dim_id))
                                S.slicing_chain(k:end) = [];
                                break
                            end
                        end
                    end
                case 'filter'
                    % filtering in dimension(s) dim_id has changed
                    k_filt = fn_find(chg_filter, [S.active_filters.obj], 'first');
                    if isempty(k_filt)
                        % filter is not active
                        if ~any(chg_filter == [S.filters.obj]), error programming, end
                        return
                    end
                    if strcmp(chg_flag, 'all')
                        % no smart update possible: filter has been
                        % completely replaced
                        S.slicing_chain(k_filt:end) = [];
                        k_start = k_filt;
                    else
                        % smart update is possible: do the smart update for
                        % this filter here                      
                        filt_k = S.active_filters(k_filt).obj;
                        % starting point data, and on which dimensions to
                        % operate
                        if k_filt == 1
                            data_start = S.data;
                        else
                            data_start = S.slicing_chain(k_filt-1);
                        end
                        % perform operation to obtain only the new part of
                        % the data that needed to be calculated
                        if isempty(ind1)
                            % no new calculation necessary: only removals
                            % and permutations
                            data_sub = [];
                        else
                            data_sub = filt_k.slicing(data_start.data, data_start.dimension_number(chg_dim_id), ind1);
                        end
                        % update the xdata object
                        chg_dim_id = filt_k.get_dim_idout(chg_dim_id); % dimension identifier in the data -> in the slice
                        head_dim = xplr.DimHeader(filt_k.headerout, chg_dim_id);                        
                        S.slicing_chain(k_filt).update_data(chg_flag, chg_dim_id, chg_ind, data_sub, head_dim)
                        % if S.slicing_chain(k_filt) is S.slice, we are done
                        if k_filt == n_active_filters
                        	return
                        end
                        k_start = k_filt + 1;
                    end
                case 'slicer'
                    % normally the invalid part of the slicing chain has
                    % already been removed, continue from there; no smart
                    % update is possible because the change is about filter
                    % addition, removal or replacement
                    k_start = length(S.slicing_chain) + 1;
            end
            % starting point
            if k_start == 1
                res = S.data;
            else
                res = S.slicing_chain(k_start - 1);
            end
            % smart updates for successive filters that can allow it
            for k = k_start:length(S.slicing_chain)
                % propagate changes in data along the slicing chain
                switch chg_flag
                    case {'remove', 'perm'}
                        % easy
                        S.slicing_chain(k).update_data(chg_flag, chg_dim_id, chg_ind, [], head_dim)
                    case {'new', 'chg', 'chg&new', 'chg&rm'}
                        % do a new slicing for a subset of the data
                        dim_idk = S.active_filters(k).dim_id;
                        filt_k = S.active_filters(k).obj;
                        % slice and update
                        dimk = res.dimension_number(dim_idk);
                        data_sub = filt_k.slicing(data_sub, dimk);
                        S.slicing_chain(k).update_data(chg_flag, chg_dim_id, chg_ind, data_sub, head_dim) % last element is the slice
                    otherwise
                        error('data_changeSmart cannot be called with flag ''%s''', chg_flag)
                end
                if k == n_active_filters
                    % the slice has been updated, we are done
                    return
                end 
                res = S.slicing_chain(k);
            end
            
            % 2) Regular slicing for the remaining part of the chain
            % separate end 
            if ~isempty(S.slicing_chain) && res == S.slice
                % res was previously the end of the chain and therefore
                % equal to the slice; but now a new filter has been
                % inserted after and we need the two xdata objects to be
                % different, otherwise the update_data_dim/chgData calls
                % applied to the slice below will also apply to this chain
                % element and this will lead to awkward results
                res = S.slice.copy();
                S.slicing_chain(end) = res;
            end
            % apply active filters one by one
            for k = length(S.slicing_chain)+1:n_active_filters
                filt_k = S.active_filters(k);
                dim_idk = filt_k.dim_id;
                obj_k = filt_k.obj;
                if ~isempty(chg_dim_id) && any(ismember(dim_idk, chg_dim_id))
                    % dimension identified by chg_dim_id in the data
                    % disappears in the slice; get the identifier of the
                    % correspond_ing dimension in the slice
                    chg_dim_id = obj_k.get_dim_idout(dim_idk);
                end
                res = obj_k.operation(res, dim_idk);
                S.slicing_chain(k) = res;
            end
            % update slice 
            switch chg_flag
                case 'global'
                    S.slice.update_data_dim('global', [], res.data, res.header)
                case 'chg_dim'
                    chg_dim_out = res.dimension_number(chg_dim_id);
                    S.slice.update_data_dim(chg_flag, chg_dim_out, res.data, res.header(chg_dim_out))
                case 'chgdata'
                    S.slice.chgData(res.data)
                otherwise
                    chg_dim_out = res.dimension_number(chg_dim_id);
                    S.slice.update_data(chg_flag, chg_dim_out, chg_ind, res.data, res.header(chg_dim_out))
            end
            % replace last chain element by slice (except if chain is
            % empty), so update_data calls to this last chain element
            % directly affect the slice (as occurs in the cases above code
            % blocks end_ing with 'return')
            if n_active_filters > 0
                S.slicing_chain(end) = S.slice;
            end
            
        end
    end
    methods
        function data_change(S,e)
            switch e.flag
                case {'global', 'chg_dim', 'all', 'new', 'chg', 'remove', 'chg&new', 'chg&rm', 'perm'}
                    % header in dimension dim_id has changed (not necessarily,
                    % but potentially also in the case 'chg'), therefore
                    % filter in that dimension is not valid anymore
                    S.slicing_chain(:) = []; % dRata has changed, all previous slicing steps became invalid
                    % check which filters remain valid
                    keep_filter = false(1, length(S.filters));
                    switch e.flag
                        case 'global'
                            % remove all filters
                        case 'chg_dim'
                            % remove filters intersecting with e.dim or
                            % whose dimension of application is not present
                            % in the data any more
                            id_x_valid = S.get_filter_index([S.data.header.dim_id]);
                            keep_filter(id_x_valid) = true;
                            id_x_dim = S.get_filter_index(e.dim);                            
                            keep_filter(id_x_dim) = false;
                        otherwise
                            % remove filters intersecting with e.dim
                            keep_filter(:) = true;
                            id_x_dim = S.get_filter_index(e.dim);                            
                            keep_filter(id_x_dim) = false;
                    end
                    % remove invalid filters: we must call viewcontrol
                    % method so that filter display is updated as well;
                    % TODO: better way of synchronizing slicer and
                    % viewcontrol's filters display
                    rmdim = [S.filters(~keep_filter).dim_id];
                    S.V.C.dimaction('rm_filter',rmdim)
                    % reslice
                    do_slice(S,'data','global')
                case 'chgdata'
                    % header was not changed, it is possible therefore to
                    % keep all existing filters
                    S.slicing_chain(:) = [];
                    do_slice(S,'data','chgdata',[],[])
                case 'name'
                    % name is not passed to the slice, so nothing to do
                otherwise
                    error 'not implemented yet'
            end
        end
        function filter_change(S,F,dim_id,e)
            if ~strcmp(e.type,'filter'), return, end
            do_slice(S,'filter',e.flag,dim_id,e.ind,F)
        end
        function reslice(S)
            S.slicing_chain(:) = [];
            do_slice(S,'data','chg_dim',1:S.nd_data,[])
        end
    end
    
end
