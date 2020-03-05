classdef filter < xplr.dataoperand
    % function F = filter(headerin,type[,label])
   
    properties (SetAccess='private')
        % input: headerin is already a property of the dataoperand mother class        
        % operation:
        type = 'selection';  % 'selection', 'indices' or 'all'
        slicefun = @nmean;   % 'nmean', 'mean', 'max', 'min', etc.
        spec
        indices = cell(1,0);
        % output: headerout is already a property of the dataoperand mother class
    end
    properties(Dependent, SetAccess='private')
        nsel
    end
    
    % Setting and updating filter
    methods
        function F = filter(headerin,type,label)
            % size and header of the input space
            if ~isa(headerin,'xplr.header'), error 'first argument must be an xplr.header object', end
            F.headerin = headerin;
            
            % operation type
            F.type = type;
            switch type
                case 'selection'
                    F.spec = xplr.selectionnd.empty(1,0);
                case {'indices' 'all'}
                    % no need for further specification
                otherwise
                    error('unknown filtering type ''%s''',type)
            end
            
            % header of the output space
            if nargin<3, label = headerin.label; end
            switch type
                case 'selection'
                    % output header is categorical, values will keep track
                    % of the selections made
                    error 'not implemented yet'
                case 'indices'
                    % output header is categorical
                    if ~isscalar(headerin)
                        % header output will be a mere enumeration (no
                        % values)
                        F.headerout = xplr.header(label,0);
                    elseif headerin.ncolumn>0
                        % categorical header: we will keep track of values
                        F.headerout = xplr.header(label,headerin.sublabels,cell(0,headerin.ncolumn));
                    elseif headerin.categorical
                        % categorical header with no values: keep track of indices
                        F.headerout = xplr.header(label,xplr.dimensionlabel('Index','numeric'),cell(0,1));
                    else
                        % measure header: keep track of values
                        F.headerout = xplr.header(label,headerin.sublabels,cell(0,1));
                    end
                case 'all'
                    % no output header because data is averaged to a single
                    % value!
                    F.headerout = xplr.header.empty(1,0);
            end
        end
        function updateSelection(F,varargin)
            % function updateSelection(F,value)
            % function updateSelection(F,'new|all',value[,'label1',headervalues1,...])
            % function updateSelection(F,'new|chg|add',ind,value[,'label1',headervalues1,...])
            % function updateSelection(F,'chg|add',ind,value)
            % function updateSelection(F,'remove|perm',ind)
            % function updateSelection(F,'reset')
            
            % check filter type
            if ~fn_ismemberstr(F.type,{'selection' 'indices'})
                error('method ''updateSelection'' is not valid for filter type ''%s''',F.type)
            end
            indicesonly = strcmp(F.type,'indices');
            
            % input
            if isscalar(varargin)
                if ischar(varargin{1})
                    if strcmp(varargin{1},'reset')
                        flag = 'all';
                        value = [];
                    else
                        error 'only flag ''reset'' can be used without arguments'
                    end
                else
                    flag = 'all';
                    value = varargin{1};
                end
                addheaderinfo = cell(2,0);
            else
                flag = varargin{1};
                varargin(1) = [];
                value = [];
                switch flag
                    case 'all'
                        value = varargin{1};
                        addheaderinfo = reshape(varargin(2:end),2,[]);
                    case {'new' 'chg' 'add' 'chg&new' 'chg&rm'}
                        if strcmp(flag,'new') && mod(nargin,2)==1
                            % flag 'new', ind not specified
                            value = varargin{1};
                            ind = F.nsel+(1:length(value));
                            addheaderinfo = reshape(varargin(2:end),2,[]);
                        else
                            [ind value] = deal(varargin{1:2});
                            addheaderinfo = reshape(varargin(3:end),2,[]);
                        end
                    case {'remove' 'perm'}
                        ind = varargin{1};
                    case 'reset'
                    otherwise
                        error('flag ''%s'' not handled by xplr.filter.updateSelection')
                end
            end
            
            
            % compute indices
            if fn_ismemberstr(flag,{'all' 'new' 'chg' 'add' 'chg&new' 'chg&rm'})
                if indicesonly
                    dataind = value;
                    if ~iscell(dataind)
                        if isscalar(dataind)
                            dataind = {dataind};
                        elseif strcmp(flag,'all')
                            % warning: the value is ambiguous, for example
                            % does value [1 2 3] mean that we want 3
                            % singleton selections (replace by {1 2 3} or
                            % a single selection of 3 indices (replace by
                            % {[1 2 3]}? We assume the first case.
                            dataind = num2cell(dataind);
                        elseif length(dataind)==length(ind)
                            dataind = num2cell(dataind);
                        else
                            error 'value should be a cell array of arrays (of indices)'
                        end
                    end
                else
                    value = value.ComputeInd([F.headerin.datalen]);
                    dataind = {value.dataind};
                end
            end
            
            % update selection
            switch flag
                case 'all'
                    % note that even if value is already equal to
                    % F.selections, we cannot just return, because
                    % addheaderinfo might bear some changes
                    ind = 1:length(dataind);
                    if ~indicesonly
                        F.selections = row(value); % in particular row(...) transform 0x0 array in 1x0 array
                    end
                    F.indices = row(dataind);
                case {'new' 'chg'}
                    if ~indicesonly, F.selections(ind) = value; end
                    F.indices(ind) = dataind;
                case 'add'
                    if ~isscalar(ind), error 'only one selection at a time can be augmented with flag ''add''', end
                    if ~indicesonly
                        F.selections(ind) = F.selections(ind).union(value); % automatic union of indices as well
                        F.indices(ind) = F.selections(ind).dataind;
                    else
                        F.indices(ind) = union(F.indices,dataind);
                    end
                    flag = 'chg'; % for the notification
                case 'chg&new'
                    if ~indicesonly, F.selections([ind{:}]) = value; end
                    F.indices([ind{:}]) = dataind;
                case 'chg&rm'
                    if ~indicesonly, F.selections(ind{1}) = value; end
                    F.indices(ind{1}) = dataind;
                    if ~indicesonly, F.selections(ind{2}) = []; end
                    F.indices(ind{2}) = [];
                case 'remove'
                    if ~indicesonly, F.selections(ind) = []; end
                    F.indices(ind) = [];
                case 'perm'
                    if ~indicesonly, F.selections = F.selections(ind); end
                    F.indices = F.indices(ind);
            end
            
            % update header of output space
            if fn_ismemberstr(flag,{'all' 'new' 'chg' 'chg&new' 'chg&rm'})
                % track header values
                switch flag
                    case 'chg&new'
                        headvalue = sliceHeader(F,[ind{:}],addheaderinfo);
                    case 'chg&rm'
                        headvalue = sliceHeader(F,ind{1},addheaderinfo);
                    otherwise
                        headvalue = sliceHeader(F,ind,addheaderinfo); % if flag is 'all', ind was set to 1:F.n
                end
                F.headerout = updateHeader(F.headerout,flag,ind,headvalue);
            else
                F.headerout = updateHeader(F.headerout,flag,ind);
            end
            
            % notification
            notify(F,'ChangedOperation',xplr.eventinfo('filter',flag,ind,value))
        end
        function setFun(F,fun)
            % convert char to function handle
            if ischar(fun)
                switch fun
                    case {'mean' 'median' 'mode' 'rms'}
                        fun = str2func(fun);
                    case 'max'
                        fun = @(x,dim)max(x,[],dim);
                    case 'min'
                        fun = @(x,dim)min(x,[],dim);
                    case 'std'
                        fun = @(x,dim)std(x,[],dim);
                    case 'var'
                        fun = @(x,dim)var(x,[],dim);
                    otherwise
                        error('unknown function for slicing ''%s''',fun)
                end
            end
            % update slicing function
            F.slicefun = fun;
            % notification
            notify(F,'ChangedOperation',xplr.eventinfo('filter','chg',1:F.nsel))
        end
    end
    
    % Get/Set Dependent
    methods
        function nsel = get.nsel(F)
            if strcmp(F.type,'all')
                nsel = 1;
            else
                nsel = length(F.indices);
            end
        end
    end
    
    % Slicing
    methods (Access='private')
        function headvalue = sliceHeader(F,ind,addheaderinfo)
            if ~fn_ismemberstr(F.type,{'selection' 'indices'})
                error('no output header for filter of type ''%s''',F.type)
            end
            
            % columns to which values have been assigned
            nind = length(ind);
            headvalue = cell(nind,F.headerout.ncolumn);
            okcolumn = false(1,F.headerout.ncolumn);
            
            % tracking of input header values
            if strcmp(F.type,'selection')
                error 'not implemented yet'
            elseif ~isscalar(F.headerin)
                % no tracking when input space has more than one dimensions
            elseif F.headerin.ncolumn>0
                % track values: call to a xplr.header method
                headvalue(:,1:F.headerin.ncolumn) = F.headerin.trackValues(F.indices(ind));
                okcolumn(1:F.headerin.ncolumn) = true;
            elseif F.headerin.categorical
                % categorical with no values: keep track of indices
                headvalue(:,1) = F.indices(ind);
                okcolumn(1) = true;
            else
                % measure header: keep track of values
                for i=1:nind
                    str = F.headerin.getItemNames(F.indices{ind(i)});
                    if isscalar(str), str = str{1}; else str = row(str); end
                    headvalue{i,1} = str;
                end
                okcolumn(1) = true;
            end
            
            % additional values set together with the selections
            if nargin<3, addheaderinfo = cell(2,0); end
            [headvalue affectedcolumns] = setAddHeaderInfo(F,headvalue,addheaderinfo);
            okcolumn(affectedcolumns) = true;
            
            % put default values in columns which were not set
            for i=find(~okcolumn)
                [headvalue{:,i}] = deal(F.headerout.sublabels(i).defaultval);
            end
        end
    end
    methods
        function slic = slicing(F,dat,dims,selsubidx)
            % slice data in dimension dim according to filter F
            % here 'dat' and 'slic' are simple Matlab arrays
            
            % input
            if nargin<4, selsubidx = 1:F.nsel; end
            nselslice = length(selsubidx);
            
            % size
            s = size(dat);
            nddata = max(max(dims),length(s));
            s(end+1:nddata) = 1;
            
            % initial reshape
            dbef = 1:min(dims)-1;
            ok = true(1,nddata); ok(dims) = false; ok(1:min(dims))=false; daft = find(ok); % faster than setdiff
            dat = fn_reshapepermute(dat,{dbef dims daft}); % this won't duplicate the array if dims are consecutive
            
            % slicing
            if nselslice==1
                switch F.type
                    case 'all'
                        slic = F.slicefun(dat,2);
                    otherwise
                        slic = F.slicefun(dat(:,F.indices{selsubidx},:),2);
                end
            else
                datatype = fn_switch(class(dat),'double','double','single');
                slic = zeros([prod(s(dbef)) nselslice prod(s(daft))],datatype);
                for i=1:nselslice
                    slic(:,i,:) = F.slicefun(dat(:,F.indices{selsubidx(i)},:),2);
                end
            end
            
            % final reshape
            sfinal = s;
            sfinal(dims(1)) = nselslice;
            sfinal(dims(2:end)) = 1;
            slic = reshape(slic,sfinal);            
        end
        function slice = operation(F,x,dims)
            % here 'data' and 'slice' are xplr.xdata objects!
            
            % check input
            checkdata(F,x,dims)
            % slice
            slic = F.slicing(x.data,dims);
            % header
            head = x.header;
            ndx = length(head);
            head = [head(1:dims(1)-1) F.headerout head(setdiff(dims(1):ndx,dims))];
            % output
            slice = xplr.xdata(slic,head);
        end
        function updateOperation(F,x,dims,slice,flag,ind)
            % here 'data' and 'slice' are xplr.xdata objects!
            
            % check input
            checkdata(F,x,dims)
            
            % slice
            switch flag
                case 'all'
                    slic = F.slicing(x.data,dims);
                case {'new' 'chg'}
                    slic = F.slicing(x.data,dims,ind);
                case {'remove' 'perm'}
                    slic = [];
                case 'chg&new'
                    slic = F.slicing(x.data,dims,[ind{:}]);
                case 'chg&rm'
                    slic = F.slicing(x.data,dims,ind{1});
                otherwise
                    error('flag ''%s'' not handled',flag)
            end
            slice.updateData(flag,dims,ind,slic,F.headerout); % this will trigger automatic notifications
        end
    end
    
    
end