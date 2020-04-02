classdef filterAndPoint < xplr.dataOperand
    % function F = filter(headerin[,label])
    % function F = filter(F,P)
    %---
    % The class filerAndPoint combines a filter (operating on a
    % n-dimensional space) and a set of n points. 
    
    properties
        F
        P = xplr.point.empty(1,0);
    end
    properties (SetAccess='protected', Transient)
        dolistenpoint = true;
    end
    properties (Dependent, SetAccess='protected', Transient)
        indices     % current selection
        index0      % current point
    end
    properties (Dependent, Transient)
        index       % current point
    end

    % Constructor and destructor
    methods
        function F = filterAndPoint(varargin)
            % init array of filters?
            if nargin==0
                return
            elseif isscalar(varargin) && isnumeric(varargin{1})
                n = varargin{1};
                if n == 0
                    F = xplr.filterAndPoint.empty(1,0);
                else
                    F(n) = xplr.filterAndPoint();
                end
                return
            end
            
            % set filter
            create_members = ~isa(varargin{1},'xplr.filter');
            if create_members
                F.F = xplr.filter(varargin{:});
            else
                F.F = varargin{1};
            end
            F.headerin = F.F.headerin;
            connectlistener(F.F,F,'ChangedOperation',@(u,e)transitNotification(F,'filter',e))
            
            % set P
            if create_members
                for i=1:F.ndin
                    F.P(i) = xplr.point(F.headerin(i));
                end
            else
                F.P = varargin{2};
                if length(F.P) ~= F.ndin, error 'number of point filters does not match number of input dimensions', end
                for i=1:F.ndin
                    if ~is_equal(F.P(i).headerin,F.headerin(i))
                        error 'headers of filter and point filter(s) do not match'
                    end
                end
            end
            for i=1:F.ndin
                connectlistener(F.P(i),F,'ChangedOperation',@(u,e)transitNotification(F,'point',e))
            end
            
            % set output header (uses filters or P depending on whether
            % filters selection is non-empty)
            setHeaderout(F)
        end
    end
    
    % Header
    methods
        function setHeaderout(F)
            F.headerout = F.F.headerout;
            if F.F.nsel==0
                % set header values
                % (tracking of input header values)
                if ~isscalar(F.headerin) 
                    disp('setting headerout for ND filter with only minimal information')
                    ncolin = 1;
                    headvalues = {fn_strcat({F.P.index},',')};
                elseif F.headerin.ncolumn>0
                    headvalues = F.headerin.values(F.P.index,:);
                    ncolin = size(headvalues,2);
                elseif F.headerin.categorical
                    ncolin = 1;
                    headvalues = {F.P.index};
                else
                    ncolin = 1;
                    headvalues = F.headerin.getItemNames(F.P.index);
                end
                % (point selection cannot define additional values: put
                % default values if such additional values were defined
                % during earlier ROI selection)
                for i=ncolin+1:F.headerout.ncolumn
                    headvalues{1,i} = F.headerout.sublabels(i).defaultval;
                end
                % replace output header
                F.headerout = updateHeader(F.headerout,'new',1,headvalues);
            end
        end
        function augmentHeader(F,varargin)
            augmentHeader(F.F,varargin{:})
            setHeaderout(F)
        end
        function [headvalue, affectedcolumns] = setAddHeaderInfo(F,~,~) %#ok<INUSD,STOUT>
            error 'method ''setAddHeaderInfo'' is not valid for filterAndPoint object; call it rather on the children filter object'
        end
    end
    
    % Point modification
    methods
        function dataind = get.indices(F)
            if F.F.nsel==0
                dataind = {fn_indices(F.szin,[F.P.index])};
            else
                dataind = F.F.indices;
            end
        end
        function x = get.index0(F)
            x = [F.P.index0];
        end
        function x = get.index(F)
            x = [F.P.index];
        end
        function set.index(F,x)
            if isequal(x,F.index0), return, end
            
            % if there are multiple points (i.e. input space is
            % multi-dimensional), it is better to generate a single event
            % after all points have been modified rather than transit one
            % event per modified point
            F.dolistenpoint = false;
            c = onCleanup(@()set(F,'dolistenpoint',true));
            
            % update point
            curind = F.index;
            for i=length(F.P), F.P(i).index = x(i); end
            
            % update header
            setHeaderout(F)
            
            % generation of a single event
            chgij = ~isequal(F.index,curind);
            notify(F,'ChangedOperation',xplr.eventinfo('point',chgij))
            if chgij && F.F.nsel==0
                notify(F,'ChangedOperation',xplr.eventinfo('filter','point'))
            end
        end
    end
    
    % Selection modification 
    methods
        function updateSelection(F,varargin)
            updateSelection(F.F,varargin{:})
            setHeaderout(F)
        end
    end
    
    % Notification
    methods
        function transitNotification(F,flag,e)
            switch flag
                case 'filter'
                    setHeaderout(F)
                    if strcmp(e.flag,'remove') && F.F.nsel==0
                        % removal of selection(s) replaces filter selection
                        % by point selection
                        e2 = xplr.eventinfo('filter','all');
                    elseif strcmp(e.flag,'new') && ismember(1,e.ind)
                        % new selection(s) replaces point selection by filter
                        % selection
                        e2 = xplr.eventinfo('filter','all');
                    else
                        % avoid very weird Matlab bug when simply passing
                        % e: make a copy of it
                        e2 = xplr.eventinfo('filter',e.flag,e.ind,e.value);
                    end
                    notify(F,'ChangedOperation',e2)
                case 'point'
                    if F.dolistenpoint % 'dolistenpoint' property is manipulated in F.set.index
                        % same as above: copy e
                        e2 = xplr.eventinfo('point',e.chgij);
                        notify(F,'ChangedOperation',e2)
                        if e.chgij
                            notify(F,'ChangedOperation',xplr.eventinfo('filter','point'))
                        end
                    end
            end
        end
    end
    
    % Slicing
    methods
        function slic = slicing(F,dat,dims,selsubidx)
            if F.F.nsel>0
                if nargin<4, selsubidx = 1:F.F.nsel; end
                slic = slicing(F.F,dat,dims,selsubidx);
            else
                if nargin>=4 && ~isequal(selsubidx,1), error 'invalid selection indices', end
                slic = slicing(F.P,dat,dims,F.ndout);
            end
        end
    end
    methods (Access='protected')
        function slic = operation_(F,dat,dims)
            % function slic = operation_(F,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            if F.F.nsel>0
                slic = F.F.operation_(dat,dims);
            else
                ndout = 1; % specify that slicing dimension should not be removed
                slic = F.P.slicing(dat,dims,ndout);
            end
        end
        function updateOperation_(F,x,dims,slice,flag,ind)
            if F.F.nsel>0
                F.F.updateOperation_(F,x,dims,slice,flag,ind)
            else
                F.P.updateOperation_(F,x,dims,slice)
            end
            %             % slice
            %             switch flag
            %                 case 'all'
            %                     if F.F.nsel>0
            %                         slic = slicing(F.F,x.data,dims);
            %                     else
            %                         if F.ndout~=1, error 'output header should be scalar', end
            %                         slic = slicing(F.P,x.data,dims,F.ndout);
            %                     end
            %                 case {'new' 'chg'}
            %                     slic = slicing(F.F,x.data,dims,ind);
            %                 case 'chg&new'
            %                     slic = slicing(F.F,x.data,dims,[ind{:}]);
            %                 case 'chg&rm'
            %                     slic = slicing(F.F,x.data,dims,ind{1});
            %                 case 'remove'
            %                     if F.F.nsel>0
            %                         slic = [];
            %                     else
            %                         flag = 'all';
            %                         if F.ndout~=1, error 'output header should be scalar', end
            %                         slic = slicing(F.P,x.data,dims,F.ndout);
            %                     end
            %                 case 'perm'
            %                     slic = [];
            %                 case 'point'
            %                     if F.ndout~=1, error 'output header should be scalar', end
            %                     flag = 'chg'; ind = 1;
            %                     slic = slicing(F.P,x.data,dims,F.ndout);
            %                 otherwise
            %                     error('flag ''%s'' not handled',flag)
            %             end
            %             slice.updateData(flag,dims,ind,slic,F.headerout); % this will trigger automatic notifications
        end
    end
    
    % Link with selection in real world coordinates
    methods
        function selection_world = operationData2Space(F)
            error 'filterAndPoint object should not be directly connected to a wordOperand object: connect rather its child point and filter'
        end
        function updateOperationData2Space(F,WO,e)
            error 'filterAndPoint object should not be directly connected to a wordOperand object: connect rather its child point and filter'
        end
        function updateOperationSpace2Data(F,~,e)
            error 'filterAndPoint object should not be directly connected to a wordOperand object: connect rather its child point and filter'
        end
    end
   
    % Copy (see xplr.dataOperand.loadfromfile)
    methods
        function copyin(F,obj)
            F.F.copyin(obj.F)
            F.P.copyin(obj.P)
        end
    end
end