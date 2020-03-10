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
    end
    properties (Access='protected')
        slicingchain = struct('res',cell(1,0),'dimdata2slice',[]); % remember the steps of slicing
        % beware that 'dim' are the dimensions of slicing IN THE ORIGINAL
        % DATA, but not any more in the slice (dimdata2slice specifies the
        % conversion)
    end
    properties (Dependent, SetAccess='private')
        nddata
        ndslice
        nactivefilt
        
    end
    
    % Constructor, destructor, basig access and get/set dependent
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
        
        
        
        %         function obj = getfilter(S,dim) % should this function return only active filter? 
        %             idx = fn_find(dim,{S.filters.dim});
        %             obj = [S.filters(idx).obj];
        %         end
    end

    % Changing filters
    methods
        function insertFilter(S,idx,dim,newfilt,doslicing)
            % function insertFilter(S,idx,dim,newfilt[,doslicing])
            if nargin<5, doslicing = true; end
            
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
            S.filters = [S.filters(1:idx-1) struct('active',true,'dim',dim,'obj',num2cell(newfilt)) S.filters(idx:end)];
            S.slicingchain(idx:end) = [];
            for i=1:nadd
                S.addListener(newfilt(i),'ChangedOperation',@(u,e)filterchange(S,dim{i},e));
            end
            % update slice
            if doslicing
                if all([newfilt.ndout]==[newfilt.ndin])
                    doslice(S,'filter','chgdim',dimadd)
                else
                    error 'not implemented yet'
                end
            end
        end
        function addFilter(S,dim,newfilt,doslicing)
            % function addFilter(S,dim,newfilt[,doslicing])
            if nargin<4, doslicing = true; end
            insertFilter(S,[],dim,newfilt,doslicing)
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
            idxinactive = idx(~active);
            idxactive = idx(active);
            % remove the filters
            S.filters(idx) = [];
            % remove invalid part of slicing chain
            % (first remove everything from the first active dim)
            if any(idxactive)
                S.slicingchain(idxactive(1):end) = []; 
            end
            % (second remove elements of the chain corresponding to removed
            % inactive filters, but without removing what follows)
            if isempty(idxactive)
                idxinactive_before_active = idxinactive;
            else
                idxinactive_before_active = idxinactive(idxinactive<idxactive(1));
            end
            if any(idxinactive_before_active)
                S.slicingchain(idxinactive_before_active) = [];
            end
            % update slice
            if doslicing && ~isempty(idxactive)
                if all([filtrm(active).ndout]==[filtrm(active).ndin])
                    activedimrm = dimrm(active);
                    doslice(S,'filter','chgdim',activedimrm)
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
            S.addListener(newfilt,'ChangedOperation',@(u,e)filterchange(S,dim,e));
            % no need for update if the filter is not active
            if ~S.filters(idx).active
                disp 'new filter replaces inactive one, so it will itself be inactive'
                return
            end
            % remove invalid part of slicing chain
            S.slicingchain(idx(1):end) = [];
            % update slice
            if doslicing
                if newfilt.ndout==prevndout
                    doslice(S,'filter','chgdim',dim)
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
            S.slicingchain(idxfirst:end) = [];
            doslice(S,'data','chgdata')
        end
        function chgFilterActive(S,idx,val)
            S.filters(idx).active = val;
            S.activateConnection(S.filters(idx).obj,val)
            S.slicingchain(idx:end) = [];
            doslice(S,'filter','chgdim',S.filters(idx).dim)
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
        function doslice(S,chgorigin,flag,dim,ind)
            % function doslice(S,chgorigin,flag,dim[,ind])
            %---
            % argument chgorigin is 'data' or 'filter', it is usefull in
            % particular for change flags such as 'new', to determine
            % whether it is the original data that has become larger in
            % dimension dim, or whether it is the filtering in dimension
            % dim which has more selection
            if nargin<4, dim = []; end
            if nargin<5, ind = []; end
            
            % smart update of the existing part of the slicing chain
            % (indices of the change)
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
            % (where to start, and get the relevant data)
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
                        datasub = subsref(S.data.data,subs);
                        % 'smart' update will not be possible for filtering in
                        % dimension dim
                        for k=1:length(S.slicingchain)
                            filterk = S.filters(k);
                            if filterk.active && any(intersect(filterk.dim,dim))
                                S.slicingchain(k:end) = [];
                                break
                            end
                        end
                    end
                case 'filter'
                    % filter in dimension dim has changed
                    kfilt = [];
                    for k=1:length(S.slicingchain)
                        filterk = S.filters(k);
                        if filterk.active && any(intersect(filterk.dim,dim))
                            if ~filterk.active
                                xplr.debuginfo('stop', ...
                                    'change in an inactive filter should not lead to slice update! please check this')
                            end
                            kfilt = k; break
                        end
                    end
                    if isempty(kfilt)
                        kstart = length(S.slicingchain)+1; 
                    elseif strcmp(flag,'all')
                        % no smart update is possible
                        S.slicingchain(kfilt:end) = [];
                        kstart = kfilt;
                    else
                        filtk = S.filters(kfilt).obj;
                        headd = filtk.headerout;
                        if isempty(ind1)
                            datasub = [];
                        else
                            if kfilt==1
                                datasub = filtk.slicing(S.data.data,dim,ind1);
                            else
                                datasub = filtk.slicing(S.slicingchain(kfilt-1).res.data,dim,ind1);
                            end
                        end
                        S.slicingchain(kfilt).res.updateData(flag,dim,ind,datasub,headd) % last element is the slice
                        if kfilt==length(S.filters), return, end % the slice has been updated, we can stop here!
                        kstart = kfilt+1;
                    end
            end
            for k = kstart:length(S.slicingchain)
                % propagate changes in data along the slicing chain
                filtk = S.filters(k);
                if ~filtk.active
                    % change was already applied, as S.slicingchain(k).res
                    % is the same as S.slicingchain(k-1).res (or S.data if
                    % k=1)
                else
                    switch flag
                        case {'remove' 'perm'}
                            % easy
                            S.slicingchain(k).res.updateData(flag,dim,ind,[],headd)
                        case {'new' 'chg' 'chg&new' 'chg&rm'}
                            % do a new slicing for a subset of the data
                            dimk = S.filters(k).dim;
                            filtk = S.filters(k).obj;
                            if isequal(dimk,dim), error 'this case is not supposed to happen (no smart update possible for filtering in the same dimension where data has changed)', end
                            datasub = filtk.slicing(datasub,dimk);
                            S.slicingchain(k).res.updateData(flag,dim,ind,datasub,headd) % last element is the slice
                        otherwise
                            error('datachangeSmart cannot be called with flag ''%s''',flag)
                    end
                end
                if k==length(S.filters), return, end % the slice has been updated, we can stop here!
            end
            
            % regular slicing for the remaining part of the chain
            % (starting point)
            if isempty(S.slicingchain)
                s = struct('res',S.data,'dimdata2slice',1:S.nddata);
            else
                s = S.slicingchain(end);
                % if this step in the chain used to be the last one but
                % will not be the last any more, its result shall no
                % longer be the slice
                if length(S.slicingchain)<length(S.filters) && s.res==S.slice
                    s.res = s.res.copy();
                    S.slicingchain(end) = s;
                end
            end
            % (apply filters one by one)
            for k = length(S.slicingchain)+1:length(S.filters)
                filtk = S.filters(k);
                dimk = filtk.dim;
                objk = filtk.obj;
                if ~filtk.active
                    % no change in data, and therefore in s
                else
                    xk = s.res;
                    s.res = objk.operation(xk,s.dimdata2slice(dimk));
                    ok = (s.dimdata2slice~=0);
                    dimbef2aft = objk.followdims(xk.nd,dimk);
                    s.dimdata2slice(ok) = dimbef2aft(s.dimdata2slice(ok));
                end
                S.slicingchain(k) = s;
            end
            
            % update slice (and replace last chain result by updated slice)
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
            if ~isempty(S.filters), S.slicingchain(end).res = S.slice; end
        end
    end
    methods
        function datachange(S,e)
            switch e.flag
                case 'global'
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    if ~isempty(S.filters)
                        rmFilter(S,1:lengh(S.filters),false) % do not update slicing here
                    end
                    doslice(S,'data','global',[],[])
                case {'chgdim' 'all' 'new' 'chg' 'remove' 'chg&new' 'chg&rm' 'perm'}
                    % header in dimension dim has changed (not necessarily,
                    % but potentially also in the case 'chg'), therefore
                    % filter in that dimension is not valid anymore
                    S.slicingchain(:) = []; % data has changed, all previous slicing steps became invalid
                    if any(ismember([S.filters.dim],e.dim))
                        rmFilterDim(S,e.dim,false) % do not update slicing here
                    end
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
        function filterchange(S,dim,e)
            flag = e.flag;
            if strcmp(flag,'point')
                flag = 'chg'; 
                ind = 1; 
            elseif isprop(e,'ind')
                ind = e.ind;
            else
                ind = [];
            end
            doslice(S,'filter',flag,dim,ind)
        end
        function reslice(S)
            S.slicingchain(:) = [];
            doslice(S,'data','chgdim',1:S.nddata)
        end
    end
    
end



