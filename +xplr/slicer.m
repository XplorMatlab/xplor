classdef slicer < xplr.graphnode
    % function S = slicer(data[,dims,filters])
    % Compute and automatically update operations applied to an initial
    % 'data', resulting in an output data called 'slice'.
    % In addition to filters, points define positions where to extract the
    % data in some dimensions in the event where there is no filter in this
    % dimension, or the filter is empty.
    
    properties 
        data
        slice
        filters = struct('active',[],'dim',cell(1,0),'obj',[]);
        % beware that 'dim' are the dimensions of slicing IN THE ORIGINAL
        % DATA, but not any more in the slice (dimdata2slice specifies the
        % conversion)
    end
    properties (SetAccess='protected')
        slicingchain = struct('res',cell(1,0),'dimdata2slice',[]); % intermediary and last steps of slicing, length is same as S.activefilters
    end
    properties (Dependent, SetAccess='private')
        nddata
        ndslice
        nactivefilt
        activefilters
    end
    
    % Constructor, destructor, basic access and get/set dependent
    methods
        function S = slicer(data,dim,filters)
            % set data
            S.data = data;
            S.addListener(data,'ChangedData',@(u,e)datachange(S,e));
            % without any filter, slice is identical data
            S.slice = data.copy();

            % set filters
            if nargin>=2 && ~isempty(filters)
                addFilter(S,dim,filters)
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
        
        
        %         function obj = getfilter(S,dim) % should this function return only active filter? 
        %             idx = fn_find(dim,{S.filters.dim});
        %             obj = [S.filters(idx).obj];
        %         end
    end

    % Changing filters
    methods
        function insertFilter(S,idx,dim,newfilt,active,doslicing)
            % function insertFilter(S,idx,dim,newfilt[active[,doslicing]])
            if nargin<5, active = true; end
            if nargin<6, doslicing = true; end
            
            % check filter
            if ~isa(newfilt,'xplr.dataoperand'), error 'filter must be a dataoperand object', end
            
            % check dimensions and headers
            if ~iscell(dim)
                if isscalar(newfilt), dim = {dim}; else dim = num2cell(dim); end
            end
            nadd = length(dim);
            dimadd = [dim{:}];
            if length(unique(dimadd))<length(dimadd)
                error 'same dimension appears in multiple filters'
            elseif any(ismember(dimadd,[S.filters.dim]))
                error 'some filter is already applied in at least one of the specified dimensions'
            elseif ~isequal([newfilt.headerin],S.data.header(dimadd))
                error 'some filter header does not match data'
            end
            % check header
            % add filters
            if isempty(idx), idx = length(S.filters)+1; end
            S.filters = [S.filters(1:idx-1) struct('active',num2cell(active),'dim',dim,'obj',num2cell(newfilt)) S.filters(idx:end)];
            % remove invalid part of the slicing chain
            nok = sum([S.filters(1:idx-1).active]);
            S.slicingchain(nok+1:end) = [];
            % install listeners
            for i=1:nadd
                S.addListener(newfilt(i),'ChangedOperation',@(u,e)filterchange(S,newfilt(i),dim{i},e));
            end
            % update slice
            if doslicing
                if all([newfilt.ndout]==[newfilt.ndin])
                    doslice(S,'slicer','chgdim',dimadd)
                else
                    error 'not implemented yet'
                end
            end
        end
        function addFilter(S,dim,newfilt,active,doslicing)
            % function addFilter(S,dim,newfilt[,active[,doslicing]])
            if nargin<4, active = true; end
            if nargin<5, doslicing = true; end
            insertFilter(S,[],dim,newfilt,active,doslicing)
        end
        function rmFilter(S,idx,doslicing)
            % function rmFilter(S,idx[,doslicing])
            if nargin<3, doslicing = true; end
            filtrm = [S.filters(idx).obj];
            S.disconnect(filtrm)
            dimrm = [S.filters(idx).dim];
            % among the dimensions for which filters will be removed, which
            % are those where the filters were really active
            active = [S.filters(idx).active];
            idxactive = idx(active);
            % remove the filters
            S.filters(idx) = [];
            % remove invalid part of slicing chain
            if ~isempty(idxactive)
                nok = sum([S.filters(1:idxactive(1)-1).active]);
                S.slicingchain(nok+1:end) = [];
            end
            % update slice
            if doslicing && ~isempty(idxactive)
                if all([filtrm(active).ndout]==[filtrm(active).ndin])
                    activedimrm = dimrm(active);
                    doslice(S,'slicer','chgdim',activedimrm)
                else
                    error 'not implemented yet'
                end
            end
        end
        function rmFilterDim(S,dim,doslicing)
            % function rmFilterDim(S,dim[,doslicing])
            if nargin<3, doslicing = true; end
            % remove filters
            idxrm = false(1,length(S.filters));
            for i=1:length(S.filters)
                idxrm(i) = any(ismember(S.filters(i).dim,dim));
            end
            if ~any(idxrm)
                fprintf('there is already no filter in dimension %i! aborting rmFilterDim\n',dim)
                return
            end
            rmFilter(S,idxrm,doslicing)
        end
        function replaceFilterDim(S,dim,newfilt,doslicing)
            % function replaceFilterDim(S,dim,newfilt[,doslicing])
            if nargin<4, doslicing = true; end
            % replace filter
            idx = fn_find(dim,{S.filters.dim});
            if isempty(idx), error 'there is no filter in the specified dimension', end
            prevndout = S.filters(idx).obj.ndout;
            S.disconnect(S.filters(idx).obj)
            S.filters(idx).obj = newfilt;
            S.addListener(newfilt,'ChangedOperation',@(u,e)filterchange(S,newfilt,dim,e));
            % no need for update if the filter is not active
            if ~S.filters(idx).active
                return
            end
            % remove invalid part of slicing chain
            nok = sum([S.filters(1:idx(1)-1).active]);
            S.slicingchain(nok+1:end) = [];
            % update slice
            if doslicing
                if newfilt.ndout==prevndout
                    doslice(S,'slicer','chgdim',dim)
                else
                    error 'not implemented yet'
                end
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
            val = logical(val);
            S.filters(idx).active = val;
            S.activateConnection(S.filters(idx).obj,val)
            nok = sum([S.filters(1:idx-1).active]);
            S.slicingchain(nok+1:end) = [];
            doslice(S,'slicer','chgdim',S.filters(idx).dim)
        end
    end

    % Get filter
    methods
        function F = getFilter(S,dim)
            idx = fn_find(dim,{S.filters.dim});
            F = [S.filters(idx).obj];
        end
    end
    
    % Slicing
    methods (Access='protected')
        function doslice(S,chgorigin,flag,dim,ind,filter)
            % function doslice(S,chgorigin,flag[,dim[,ind[,filter]]])
            %---
            % Apply filters successively to get a slice out of the data.
            % "Smart updates" apply when the flag (for example 'new')
            % indicates that only a part of the data will need to be
            % re-calculated (for example with 'new', all the existing slice
            % will be preserved, but some additional part will be added
            % corresponding to the original increase in data or filter).

            if nargin<4, dim = []; end
            if nargin<5, ind = []; end
            if nargin<6, filter = []; end
            
            nactivefilters = sum([S.filters.active]);
            
            % 1) Smart update of the existing part of the slicing chain
            % indices of the change
            switch flag
                case {'new' 'chg'}
                    ind1 = ind;
                case 'chg&new'
                    ind1 = [ind{:}];
                case 'chg&rm'
                    ind1 = ind{1};
                otherwise
                    ind1 = [];
            end
            % where to start, and get the relevant sub-part of the data
            % needed for smart update calculations
            switch chgorigin
                case 'data'
                    kstart = 1;
                    if fn_ismemberstr(flag,{'all' 'chgdata' 'global' 'chgdim' 'insertdim' 'rmdim' 'permdim'})
                        % smart update is not possible
                        S.slicingchain(:) = [];
                    else
                        % original data has changed in dimension dim
                        headd = S.data.header(dim);
                        subs = substruct('()',repmat({':'},1,S.data.nd));
                        subs.subs{dim} = ind1;
                        % part of the data that will be used for new
                        % calculations, while preserving some of the
                        % existing calculations
                        datasub = subsref(S.data.data,subs); 
                        % the 'smart' update however will not be possible
                        % once there will be another filter applied in the
                        % same dimension dim
                        for k=1:length(S.slicingchain)
                            filterk = S.activefilters(k);
                            if filterk.active && any(intersect(filterk.dim,dim))
                                S.slicingchain(k:end) = [];
                                break
                            end
                        end
                    end
                case 'filter'
                    % filtering in dimension(s) dim has changed
                    kfilt = fn_find([S.activefilters.obj],filter,'first');
                    if isempty(kfilt)
                        % filter is not active
                        if ~any(filter == [S.filters.obj]), error programming, end
                        return
                    end
                    if strcmp(flag,'all')
                        % no smart update possible: filter has been
                        % completely replaced
                        S.slicingchain(kfilt:end) = [];
                        kstart = kfilt;
                    else
                        % smart update is possible: do the smart update for
                        % this filter here                      
                        filtk = S.activefilters(kfilt).obj;
                        headd = filtk.headerout;
                        if isempty(ind1)
                            % no new calculation necessary: only removals
                            % and permutations
                            datasub = [];
                        else
                            if kfilt==1
                                datasub = filtk.slicing(S.data.data,dim,ind1);
                            else
                                datasub = filtk.slicing(S.slicingchain(kfilt-1).res.data,dim,ind1);
                            end
                        end
                        S.slicingchain(kfilt).res.updateData(flag,dim,ind,datasub,headd)
                        if kfilt == nactivefilters
                            % S.slicingchain(kfilt) is S.slice, we are done
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
            % smart updates for successive filters that can allow it
            for k = kstart:length(S.slicingchain)
                % propagate changes in data along the slicing chain
                switch flag
                    case {'remove' 'perm'}
                        % easy
                        S.slicingchain(k).res.updateData(flag,dim,ind,[],headd)
                    case {'new' 'chg' 'chg&new' 'chg&rm'}
                        % do a new slicing for a subset of the data
                        dimk = S.activefilters(k).dim;
                        filtk = S.activefilters(k).obj;
                        if isequal(dimk,dim), error 'this case is not supposed to happen (no smart update possible for filtering in the same dimension where data has changed)', end
                        datasub = filtk.slicing(datasub,dimk);
                        S.slicingchain(k).res.updateData(flag,dim,ind,datasub,headd) % last element is the slice
                    otherwise
                        error('datachangeSmart cannot be called with flag ''%s''',flag)
                end
                if k == nactivefilters
                    % the slice has been updated, we are done
                    return
                end 
            end
            
            % 2) Regular slicing for the remaining part of the chain
            % starting point
            if isempty(S.slicingchain)
                s = struct('res',S.data,'dimdata2slice',1:S.nddata);
            else
                s = S.slicingchain(end);
                if s.res==S.slice 
                    % s was previously the end of the chain and therefore
                    % equal to the slice; but now a new filter has been
                    % inserted after and we need the two xdata object to be
                    % different otherwise the updateDataDim/chgData calls
                    % below will lead to awkward results
                    s.res = S.slice.copy();
                    S.slicingchain(end) = s;
                end
            end
            % apply active filters one by one
            for k = length(S.slicingchain)+1:nactivefilters
                filtk = S.activefilters(k);
                dimk = filtk.dim;
                objk = filtk.obj;
                resprev = s.res;
                s.res = objk.operation(resprev,s.dimdata2slice(dimk));
                ok = (s.dimdata2slice~=0);
                dimbef2aft = objk.followdims(resprev.nd,dimk);
                s.dimdata2slice(ok) = dimbef2aft(s.dimdata2slice(ok));
                S.slicingchain(k) = s;
            end
            % update slice 
            switch flag
                case 'global'
                    S.slice.updateDataDim(flag,[],s.res.data,s.res.header)
                case {'chgdim' 'insertdim'}
                    S.slice.updateDataDim(flag,dim,s.res.data,s.res.header(dim))
                case {'rmdim' 'permdim'}
                    S.slice.updateDataDim(flag,dim,s.res.data)
                case 'chgdata'
                    S.slice.chgData(s.res.data)
                otherwise
                    S.slice.updateData(flag,dim,ind,s.res.data,s.res.header(dim))
            end
            % replace last chain element by slice (except if chain is
            % empty), so updateData calls to this last chain element
            % directly affect the slice (as occurs in the cases above code
            % blocks ending with 'return')
            if nactivefilters > 0
                S.slicingchain(end).res = S.slice;
            end
            
        end
    end
    methods
        function datachange(S,e)
            switch e.flag
                case 'global'
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    doslice(S,'data','global',[],[])
                case {'chgdim' 'all' 'new' 'chg' 'remove' 'chg&new' 'chg&rm' 'perm'}
                    % header in dimension dim has changed (not necessarily,
                    % but potentially also in the case 'chg'), therefore
                    % filter in that dimension is not valid anymore
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    doslice(S,'data','chgdim',e.dim,[])
                case 'chgdata'
                    % header was not changed, it is possible therefore to
                    % keep all existing filters
                    S.slicingchain(:) = [];
                    doslice(S,'data','chgdata',[],[])
                otherwise
                    error 'not implemented yet'
            end
        end
        function filterchange(S,F,dim,e)
            flag = e.flag;
            if strcmp(flag,'point')
                flag = 'chg'; 
                ind = 1; 
            else
                ind = e.ind;
            end
            doslice(S,'filter',flag,dim,ind,F)
        end
        function reslice(S)
            S.slicingchain(:) = [];
            doslice(S,'data','chgdim',1:S.nddata,[])
        end
    end
    
end



