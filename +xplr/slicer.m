classdef slicer < xplr.graphnode
    % function S = slicer(data[,dims,filters])
    % Compute and automatically update operations applied to an initial
    % 'data', resulting in an output data called 'slice'.
    % In addition to filters, points define positions where to extract the
    % data in some dimensions in the event where there is no filter in this
    % dimension, or the filter is empty.
    
    % In general, data dimensions are identified by their unique identifier
    % (dimID) rather than by their number, but some methods accept both
    % identification methods (dimension numbers can easily be distinguished
    % from dimension identifiers because the former are >=1, the latter are
    % <1 ; conversion of either number, identifier or label to identifier
    % is performed by the xplr.xdata.dimensionID method).
    
    properties 
        V
        data
        slice
        filters = struct('active',[],'dimID',cell(1,0),'obj',[]);
    end
    properties (SetAccess='protected')
        slicingchain = xplr.xdata.empty(1,0); % intermediary and last steps of slicing, length is same as S.activefilters
        pendingrmfilter = false; % remember when some dimensions were removed but reslice did not occur yet
    end
    properties (Dependent, SetAccess='private')
        nddata
        ndslice
        nactivefilt
        activefilters
    end
    
    % Constructor, destructor, basic access and get/set dependent
    methods
        function S = slicer(V,data,dimID,filters)
            % link to parent view
            S.V = V;
            xplr.debuginfo('TODO','can a slicer exist without being aware of its possessing view?')
            % set data
            S.data = data;
            S.addListener(data,'ChangedData',@(u,e)datachange(S,e));
            % without any filter, slice is identical data
            S.slice = data.copy();

            % set filters
            if nargin>=3 && ~isempty(filters)
                addFilter(S,dimID,filters)
            end
        end
        function n = get.nddata(S)
            n = S.data.nd;
        end
        function n = get.ndslice(S)
            n = S.slice.nd;
        end
        function n = get.nactivefilt(S)
            n = sum([S.filters.active]);
        end
        function filters = get.activefilters(S)
            filters = S.filters([S.filters.active]);
        end
    end

    % Changing filters
    methods
        function insertFilter(S,idx,dimID,newfilt,active)
            % function insertFilter(S,idx,dimID,newfilt[active])
            dimID = S.data.dimensionID(dimID);
            if nargin<5, active = true; end
            
            % check filter
            if ~isa(newfilt,'xplr.dataOperand'), error 'filter must be a dataOperand object', end
            
            % check dimensions and headers
            if ~iscell(dimID)
                if isscalar(newfilt), dimID = {dimID}; else dimID = num2cell(dimID); end
            end
            nadd = length(dimID);
            dimIDadd = [dimID{:}];
            if length(unique(dimIDadd))<length(dimIDadd)
                error 'same dimension appears in multiple filters'
            elseif any(ismember(dimIDadd,[S.filters.dimID]))
                error 'some filter is already applied in at least one of the specified dimensions'
            elseif ~isequal([newfilt.headerin],S.data.headerByID(dimIDadd))
                error 'some filter header does not match data'
            end
            % check header
            % add filters
            if isempty(idx), idx = length(S.filters)+1; end
            S.filters = [S.filters(1:idx-1) struct('active',num2cell(active),'dimID',dimID,'obj',num2cell(newfilt)) S.filters(idx:end)];
            % remove invalid part of the slicing chain
            nok = sum([S.filters(1:idx-1).active]);
            S.slicingchain(nok+1:end) = [];
            % install listeners
            for i=1:nadd
                S.addListener(newfilt(i),'ChangedOperation',@(u,e)filterchange(S,newfilt(i),dimID{i},e));
            end
            % update slice
            if ~S.pendingrmfilter && isscalar(newfilt) && newfilt.ndout==newfilt.ndin
                doslice(S,'slicer','chgdim',dimIDadd)
            else
                doslice(S,'slicer','global')
                S.pendingrmfilter = false;
            end
        end
        function addFilter(S,dimID,newfilt,active)
            % function addFilter(S,dimID,newfilt[,active])
            if nargin<4, active = true; end
            insertFilter(S,[],dimID,newfilt,active)
        end
        function rmFilter(S,idx,doslicing)
            % function rmFilter(S,idx[,doslicing])
            if nargin<3, doslicing = true; end
            filtrm = [S.filters(idx).obj];
            S.disconnect(filtrm)
            % among the dimensions for which filters will be removed, which
            % are those where the filters were really active
            active = [S.filters(idx).active];
            idxactive = idx(active);
            chgdimID = [S.filters(idxactive).dimID];
            % remove invalid part of slicing chain
            if ~isempty(idxactive)
                nok = sum([S.filters(1:idxactive(1)-1).active]);
                S.slicingchain(nok+1:end) = [];
            end
            % remove the filters
            S.filters(idx) = [];
            % update slice
            if isempty(idxactive), return, end
            if doslicing
                if all(isvalid(filtrm)) && all([filtrm(active).ndout]==[filtrm(active).ndin])
                    doslice(S,'slicer','chgdim',chgdimID)
                else
                    doslice(S,'slicer','global')
                end
            else
                % remember that reslice did not occur after removing filter
                S.pendingrmfilter = true;
            end
        end
        function rmFilterDim(S,dimID,doslicing)
            % function rmFilterDim(S,dimID[,doslicing])
            dimID = S.data.dimensionID(dimID);
            if nargin<3, doslicing = true; end
            % remove filters
            idxrm = S.getFilterIndex(dimID);
            if ~any(idxrm)
                fprintf('there is already no filter in dimension %i! aborting rmFilterDim\n',dimID)
                return
            end
            rmFilter(S,idxrm,doslicing)
        end
        function replaceFilter(S,idxs,dimIDs,newfilt)
            % function replaceFilter(S,idxs,dimIDs,newfilt)
            
            % there can be several filters: dimIDs must be a cell array
            dimIDs = S.data.dimensionID(dimIDs);
            if ~iscell(dimIDs)
                if ~isscalar(newfilt), dimIDs = num2cell(dimIDs); else, dimIDs = {dimIDs}; end
            end    
            if any(diff([length(idxs), length(dimIDs), length(newfilt)])), error 'length mismatch', end

            % replace filter(s)
            ndout_changed = false;
            chgdim = zeros(1,0);
            for i = 1:length(dimIDs)
                dimID = dimIDs{i};
                F = newfilt(i);
                idx = idxs(i);
                % replace filter
                prevndout = S.filters(idx).obj.ndout;
                S.disconnect(S.filters(idx).obj)
                S.filters(idx).obj = F;
                S.filters(idx).dimID = dimID;
                S.addListener(F,'ChangedOperation',@(u,e)filterchange(S,F,dimID,e));
                % no need for update if the filter is not active
                if S.filters(idx).active
                    if F.ndout~=prevndout, ndout_changed = true; end
                    chgdim = [chgdim dimID]; %#ok<AGROW>
                else
                    S.activateConnection(F,false)
                end
                % remove invalid part of slicing chain
                nok = sum([S.filters(1:idx(1)-1).active]);
                S.slicingchain(nok+1:end) = [];
            end
            
            % update slice
            if ndout_changed
                doslice(S,'slicer','global')
            else
                doslice(S,'slicer','chgdim',chgdim)
            end
        end
        function replaceFilterDim(S,dimIDs,newfilt)
            % function replaceFilterDim(S,dims,newfilt)
            
            % there can be several filters: dimIDs must be a cell array
            dimIDs = S.data.dimensionID(dimIDs);
            if ~iscell(dimIDs)
                if ~isscalar(newfilt), dimIDs = num2cell(dimIDs); else, dimIDs = {dimIDs}; end
            end               

            % replace filter(s)
            ndout_changed = false;
            chgdim = zeros(1,0);
            for i = 1:length(dimIDs)
                dimID = dimIDs{i};
                F = newfilt(i);
                % replace filter
                idx = fn_find(dimID,{S.filters.dimID});
                if isempty(idx), error 'there is no filter in the specified dimension', end
                prevndout = S.filters(idx).obj.ndout;
                S.disconnect(S.filters(idx).obj)
                S.filters(idx).obj = F;
                S.addListener(F,'ChangedOperation',@(u,e)filterchange(S,F,dimID,e));
                % no need for update if the filter is not active
                if S.filters(idx).active
                    if F.ndout~=prevndout, ndout_changed = true; end
                    chgdim = [chgdim dimID]; %#ok<AGROW>
                else
                    S.activateConnection(F,false)
                end
                % remove invalid part of slicing chain
                nok = sum([S.filters(1:idx(1)-1).active]);
                S.slicingchain(nok+1:end) = [];
            end
            
            % update slice
            if ndout_changed
                doslice(S,'slicer','global')
            else
                doslice(S,'slicer','chgdim',chgdim)
            end
        end
        function permFilters(S,perm)
            % function permFilters(S,perm)
            % modify the order in which filters are applied
            
            % check
            if ~isequal(sort(perm),1:length(S.filters))
                error 'argument is not a valid permutation'
            end
            idxfirst = find(perm~=1:length(S.filters),1,'first');
            if isempty(idxfirst)
                % no permutation at all!
                return
            end
            
            % permutation
            S.filters = S.filters(perm);
            
            % update slice: note that all output header will remain the
            % same, so best is to signal the change as a change in the data
            nok = sum([S.filters(1:idxfirst-1).active]);
            S.slicingchain(nok+1:end) = [];
            doslice(S,'slicer','chgdata')
        end
        function chgFilterActive(S,idx,val)
            % function chgFilterActive(S,idx,val)
            val = boolean(val);
            if all([S.filters(idx).active] == val), return, end
            for i = idx
                S.filters(i).active = val;
                S.activateConnection(S.filters(i).obj,val)
            end
            nok = sum([S.filters(1:min(idx)-1).active]);
            S.slicingchain(nok+1:end) = [];
            doslice(S,'slicer','global')
        end
        function applyPending(S)
            if S.pendingrmfilter
                S.doslice('slicer','global')
                S.pendingrmfilter = false;
            end
        end
    end

    % Get filter
    methods
        function F = getFilterByDim(S,dimID)
            dimID = S.data.dimensionID(dimID);
            idx = fn_find(dimID,{S.filters.dimID});
            F = [S.filters(idx).obj];
        end
        function idx = getFilterIndex(S,F_or_dim)
            % function idx = getFilterIndex(S,[F|dim|dimID])
            if isnumeric(F_or_dim)
                dimID = S.data.dimensionID(F_or_dim);
                idx = false(1,length(S.filters));
                for i=1:length(S.filters)
                    idx(i) = any(ismember(S.filters(i).dimID,dimID));
                end
                idx = find(idx);
            else
                F = F_or_dim;
                if ~isscalar(F), error 'input filter must be scalar', end
                idx = fn_find(F,[S.filters.obj],'first');
            end
        end
    end
    
    % Slicing
    methods (Access='protected')
        function doslice(S,chgorigin,chgflag,chgdim,chgind,chgfilter)
            % function doslice(S,chgorigin,flag[,chgdim[,ind[,filter]]])
            %---
            % Apply filters successively to get a slice out of the data.
            % "Smart updates" apply when the flag (for example 'new')
            % indicates that only a part of the data will need to be
            % re-calculated (for example with 'new', all the existing slice
            % will be preserved, but some additional part will be added
            % corresponding to the original increase in data or filter).

            if nargin<4, chgdim = []; end
            if nargin<5, chgind = []; end
            if nargin<6, chgfilter = []; end
            
            nactivefilters = sum([S.filters.active]);
            
            % 1) Smart update of the existing part of the slicing chain
            % dimension number
            [chgdim, chgdimID] = S.data.dimensionNumberAndID(chgdim);
            % indices of the change
            switch chgflag
                case {'new' 'chg'}
                    ind1 = chgind;
                case 'chg&new'
                    ind1 = [chgind{:}];
                case 'chg&rm'
                    ind1 = chgind{1};
                otherwise
                    ind1 = [];
            end
            % where to start, and get the relevant sub-part of the data
            % needed for smart update calculations
            switch chgorigin
                case 'data'
                    kstart = 1;
                    if fn_ismemberstr(chgflag,{'all' 'chgdata' 'global' 'chgdim'})
                        % smart update is not possible
                        S.slicingchain(:) = [];
                    else
                        if ~isscalar(chgdim), error 'programming', end
                        % original data has changed in dimension dimID
                        headdim = S.data.header(chgdim);
                        subs = substruct('()',repmat({':'},1,S.data.nd));
                        subs.subs{chgdim} = ind1;
                        % part of the data that will be used for new
                        % calculations, while preserving some of the
                        % existing calculations
                        datasub = subsref(S.data.data,subs); 
                        % the 'smart' update however will not be possible
                        % once there will be another filter applied in the
                        % same dimension dim
                        for k=1:length(S.slicingchain)
                            filterk = S.activefilters(k);
                            if filterk.active && any(intersect(filterk.dimID,chgdimID))
                                S.slicingchain(k:end) = [];
                                break
                            end
                        end
                    end
                case 'filter'
                    % filtering in dimension(s) dimID has changed
                    kfilt = fn_find(chgfilter,[S.activefilters.obj],'first');
                    if isempty(kfilt)
                        % filter is not active
                        if ~any(chgfilter == [S.filters.obj]), error programming, end
                        return
                    end
                    if strcmp(chgflag,'all')
                        % no smart update possible: filter has been
                        % completely replaced
                        S.slicingchain(kfilt:end) = [];
                        kstart = kfilt;
                    else
                        % smart update is possible: do the smart update for
                        % this filter here                      
                        filtk = S.activefilters(kfilt).obj;
                        % starting point data, and on which dimensions to
                        % operate
                        if kfilt == 1
                            datastart = S.data;
                        else
                            datastart = S.slicingchain(kfilt-1);
                        end
                        % perform operation to obtain only the new part of
                        % the data that needed to be calculated
                        if isempty(ind1)
                            % no new calculation necessary: only removals
                            % and permutations
                            datasub = [];
                        else
                            datasub = filtk.slicing(datastart.data,datastart.dimensionNumber(chgdimID),ind1);
                        end
                        % update the xdata object
                        chgdimID = filtk.getdimIDout(chgdimID); % dimension identifier in the data -> in the slice
                        headdim = xplr.dimheader(filtk.headerout,chgdimID);                        
                        S.slicingchain(kfilt).updateData(chgflag,chgdimID,chgind,datasub,headdim)
                        % if S.slicingchain(kfilt) is S.slice, we are done
                        if kfilt == nactivefilters
                        	return
                        end
                        kstart = kfilt+1;
                    end
                case 'slicer'
                    % normally the invalid part of the slicing chain has
                    % already been removed, continue from there; no smart
                    % update is possible because the change is about filter
                    % addition, removal or replacement
                    kstart = length(S.slicingchain)+1;
            end
            % starting point
            if kstart == 1
                res = S.data;
            else
                res = S.slicingchain(kstart-1);
            end
            % smart updates for successive filters that can allow it
            for k = kstart:length(S.slicingchain)
                % propagate changes in data along the slicing chain
                switch chgflag
                    case {'remove' 'perm'}
                        % easy
                        S.slicingchain(k).updateData(chgflag,chgdimID,chgind,[],headdim)
                    case {'new' 'chg' 'chg&new' 'chg&rm'}
                        % do a new slicing for a subset of the data
                        dimIDk = S.activefilters(k).dimID;
                        filtk = S.activefilters(k).obj;
                        % slice and update
                        dimk = res.dimensionNumber(dimIDk);
                        datasub = filtk.slicing(datasub,dimk);
                        S.slicingchain(k).updateData(chgflag,chgdimID,chgind,datasub,headdim) % last element is the slice
                    otherwise
                        error('datachangeSmart cannot be called with flag ''%s''',chgflag)
                end
                if k == nactivefilters
                    % the slice has been updated, we are done
                    return
                end 
                res = S.slicingchain(k);
            end
            
            % 2) Regular slicing for the remaining part of the chain
            % separate end 
            if ~isempty(S.slicingchain) && res==S.slice
                % res was previously the end of the chain and therefore
                % equal to the slice; but now a new filter has been
                % inserted after and we need the two xdata objects to be
                % different, otherwise the updateDataDim/chgData calls
                % applied to the slice below will also apply to this chain
                % element and this will lead to awkward results
                res = S.slice.copy();
                S.slicingchain(end) = res;
            end
            % apply active filters one by one
            for k = length(S.slicingchain)+1:nactivefilters
                filtk = S.activefilters(k);
                dimIDk = filtk.dimID;
                objk = filtk.obj;
                if ~isempty(chgdimID) && any(ismember(dimIDk,chgdimID))
                    % dimension identified by chgdimID in the data
                    % disappears in the slice; get the identifier of the
                    % corresponding dimension in the slice
                    chgdimID = objk.getdimIDout(dimIDk);
                end
                res = objk.operation(res,dimIDk);
                S.slicingchain(k) = res;
            end
            % update slice 
            switch chgflag
                case 'global'
                    S.slice.updateDataDim('global',[],res.data,res.header)
                case 'chgdim'
                    chgdimout = res.dimensionNumber(chgdimID);
                    S.slice.updateDataDim(chgflag,chgdimout,res.data,res.header(chgdimout))
                case 'chgdata'
                    S.slice.chgData(res.data)
                otherwise
                    chgdimout = res.dimensionNumber(chgdimID);
                    S.slice.updateData(chgflag,chgdimout,chgind,res.data,res.header(chgdimout))
            end
            % replace last chain element by slice (except if chain is
            % empty), so updateData calls to this last chain element
            % directly affect the slice (as occurs in the cases above code
            % blocks ending with 'return')
            if nactivefilters > 0
                S.slicingchain(end) = S.slice;
            end
            
        end
    end
    methods
        function datachange(S,e)
            switch e.flag
                case {'global' 'chgdim' 'all' 'new' 'chg' 'remove' 'chg&new' 'chg&rm' 'perm'}
                    % header in dimension dimID has changed (not necessarily,
                    % but potentially also in the case 'chg'), therefore
                    % filter in that dimension is not valid anymore
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    % check which filters remain valid
                    keepfilter = false(1,length(S.filters));
                    switch e.flag
                        case 'global'
                            % remove all filters
                        case 'chgdim'
                            % remove filters intersecting with e.dim or
                            % whose dimension of application is not present
                            % in the data any more
                            idxvalid = S.getFilterIndex([S.data.header.dimID]);
                            keepfilter(idxvalid) = true;
                            idxdim = S.getFilterIndex(e.dim);                            
                            keepfilter(idxdim) = false;
                        otherwise
                            % remove filters intersecting with e.dim
                            keepfilter(:) = true;
                            idxdim = S.getFilterIndex(e.dim);                            
                            keepfilter(idxdim) = false;
                    end
                    % remove invalid filters: we must call viewcontrol
                    % method so that filter display is updated as well;
                    % TODO: better way of synchronizing slicer and
                    % viewcontrol's filters display
                    rmdim = [S.filters(~keepfilter).dimID];
                    S.V.C.dimaction('rmfilter',rmdim)
                    % reslice
                    doslice(S,'data','global')
                case 'chgdata'
                    % header was not changed, it is possible therefore to
                    % keep all existing filters
                    S.slicingchain(:) = [];
                    doslice(S,'data','chgdata',[],[])
                case 'name'
                    % name is not passed to the slice, so nothing to do
                otherwise
                    error 'not implemented yet'
            end
        end
        function filterchange(S,F,dimID,e)
            if ~strcmp(e.type,'filter'), return, end
            doslice(S,'filter',e.flag,dimID,e.ind,F)
        end
        function reslice(S)
            S.slicingchain(:) = [];
            doslice(S,'data','chgdim',1:S.nddata,[])
        end
    end
    
end



