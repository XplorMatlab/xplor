classdef filter < xplr.dataOperand
    % function F = filter(headerin[,label])
   
    properties (SetAccess='protected')
        % input: headerin is already a property of the dataOperand mother class
        % operation:
        selection = xplr.selectionnd.empty(1,0);
        slicefun = @nmean;   % 'nmean', 'mean', 'max', 'min', etc.
        % output: headerout is already a property of the dataOperand mother class
    end
    properties(Dependent, SetAccess='protected', Transient)
        nsel
        indices
    end
    
    % Setting and updating filter
%     methods (Static)
%         function s = loadobj(obj)
%             s = obj;
%         end
%     end
    methods
        function F = filter(headerin,label)
            % size and header of the input space
            if ~isa(headerin,'xplr.header'), error 'first argument must be an xplr.header object', end
            F.headerin = headerin;
            
            % header of the output space
            if nargin<2
                if isscalar(headerin)
                    label = headerin.label;
                else
                    label = fn_strcat({headerin.label},'(',',',')');
                end
            end

            % output header is categorical
            if ~isscalar(headerin)
                % header output will be a mere enumeration (no values)
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
        end
        function updateSelection(F,varargin)
            % function updateSelection(F,value)
            % function updateSelection(F,'new|all',value[,'label1',headervalues1,...])
            % function updateSelection(F,'new|all|chg|add',ind,value[,'label1',headervalues1,...])
            % function updateSelection(F,'chg|add',ind,value)
            % function updateSelection(F,'remove|perm',ind)
            % function updateSelection(F,'reset')
            
            % input
            if isscalar(varargin)
                if ischar(varargin{1})
                    if strcmp(varargin{1},'reset')
                        flag = 'all';
                        value = xplr.selectionnd.empty(1,0);
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
                switch flag
                    case {'all' 'new' 'chg' 'add' 'chg&new' 'chg&rm'}
                        if mod(nargin,2)==1
                            % ind not specified
                            value = varargin{1};
                            switch flag
                                case 'all'
                                    ind = 1:length(value);
                                case 'new'
                                    ind = F.nsel+(1:length(value));
                                otherwise
                                    error 'incorrect number of arguments'
                            end
                            addheaderinfo = reshape(varargin(2:end),2,[]);
                        else
                            % ind specified
                            [ind, value] = deal(varargin{1:2});
                            addheaderinfo = reshape(varargin(3:end),2,[]);
                        end
                    case {'remove' 'perm'}
                        ind = varargin{1};
                        value = []; % should not be used
                    otherwise
                        error('flag ''%s'' not handled by xplr.filter.updateSelection')
                end
            end
            
            % compute indices
            if fn_ismemberstr(flag,{'all' 'new' 'chg' 'add' 'chg&new' 'chg&rm'})
                value = value.ComputeInd([F.headerin.n]);
                % filter definition when data is categorical can only be
                % indices based (but not shape based)
                if any([F.headerin.categorical])
                    value = value.convert('indices',[F.headerin.n]);
                end
            end
            
            % update selection
            switch flag
                case 'all'
                    % note that even if value is already equal to
                    % F.selection, we cannot just return, because
                    % addheaderinfo might bear some changes
                    ind = 1:length(value);
                    F.selection = row(value); % in particular row(...) transform 0x0 array in 1x0 array
                case {'new' 'chg'}
                    F.selection(ind) = value;
                case 'add'
                    if ~isscalar(ind), error 'only one selection at a time can be augmented with flag ''add''', end
                    F.selection(ind) = F.selection(ind).union(value); % automatic union of indices as well
                    flag = 'chg'; % for the notification
                case 'chg&new'
                    F.selection([ind{:}]) = value;
                case 'chg&rm'
                    F.selection(ind{1}) = value;
                    F.selection(ind{2}) = [];
                case 'remove'
                    F.selection(ind) = [];
                case 'perm'
                    F.selection = F.selection(ind);
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
            e = xplr.eventinfo('filter',flag,ind,value);
            notify(F,'ChangedOperation',e)
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
        function copyin(F,obj)
            % do not call updateing methods because there might be
            % additional information in headerout; change manually the
            % needed properties and raise event
            F.slicefun = obj.slicefun;
            F.selection = obj.selection;
            F.headerout = obj.headerout;
            e = xplr.eventinfo('filter','all',1:length(F.selection),F.selection);
            notify(F,'ChangedOperation',e)
        end
    end
    
    % Get/Set Dependent
    methods
        function nsel = get.nsel(F)
            nsel = length(F.selection);
        end
        function indices = get.indices(F)
            % get data indices from selections
            indices = cell(1,F.nsel);
            for i = 1:F.nsel
                indices{i} = F.selection(i).ComputeInd([F.headerin.n]).dataind;
            end
        end
    end
    
    % Slicing
    methods (Access='private')
        function headvalue = sliceHeader(F,ind,addheaderinfo)
            % columns to which values have been assigned
            nind = length(ind);
            headvalue = cell(nind,F.headerout.ncolumn);
            okcolumn = false(1,F.headerout.ncolumn);
            
            % tracking of input header values
            if ~isscalar(F.headerin)
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
            % function slic = slicing(F,dat,dims,selsubidx)
            %---
            % slice data in dimension dim according to filter F
            % dat and slic are simple Matlab arrays
            
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
                slic = F.slicefun(dat(:,F.indices{selsubidx},:),2);
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
            sfinal(dims(2:end)) = [];
            slic = reshape(slic,sfinal);            
        end
    end
    methods (Access='protected')
        function slic = operation_(F,dat,dims)
            % function slic = operation_(F,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            slic = F.slicing(dat,dims);
        end
        function updateOperation_(F,x,dims,slice,flag,ind)
            % function updateOperation_(F,x,dims,slice,flag,ind)
            %---
            % x and slice are xplr.xdata objects
            
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
    
    % Link with selection in real world coordinates
    methods
        function selection_world = operationData2Space(F)
            aff = xplr.affinitynd([F.headerin.scale],[F.headerin.start]-[F.headerin.scale]);
            selection_world = aff.move_selection(F.selection);
        end
        function updateOperationData2Space(F,WO,e)
            [flag, ind, value] = deal(e.flag, e.ind, e.value);
            if ~isempty(value)
                aff = xplr.affinitynd([F.headerin.scale],[F.headerin.start]-[F.headerin.scale]);
                value = aff.move_selection(value);
            end
            switch flag
                case 'all'
                    WO.operation = value;
                case {'new' 'chg'}
                    WO.operation(ind) = value;
                case 'remove'
                    WO.operation(ind) = [];
                case 'permute'
                    WO.operation = WO.operation(ind);
                case 'chg&new'
                    WO.operation([ind{:}]) = value;
                case 'chg&rm'
                    WO.operation(ind{1}) = value;
                    WO.operation(ind{2}) = [];
                otherwise
                    error('flag ''%s'' not handled',flag)
            end
            notify(WO,'ChangedOperation',xplr.eventinfo('filter',flag,ind,value))
        end
        function updateOperationSpace2Data(F,selection_world,e)
            if nargin>=3
                [flag, ind, value] = deal(e.flag, e.ind, e.value);
            else
                flag = 'all';
                value = selection_world;
                ind = 1:length(value);
            end
            if ~isempty(value)
                aff = xplr.affinitynd(1./[F.headerin.scale],1-[F.headerin.start]./[F.headerin.scale]);
                value = value.applyaffinity(aff);
            end
            F.updateSelection(flag,ind,value)
        end
    end
    
    
end