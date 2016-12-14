classdef filterAndPoint < xplr.dataoperand
    % The class filerAndPoint combines a filter (operating on a
    % n-dimensional space) and a set of n points. 
    
    properties
        F
        P = xplr.point.empty(1,0);
    end
    properties (SetAccess='private')
        listeners = struct('filter',[],'point',event.listener.empty(1,0),'pointselection',event.listener.empty(1,0));
        dolistenpoint
    end
    properties (Dependent, SetAccess='private')
        indices     % current selection
        index0      % current point
    end
    properties (Dependent)
        index       % current point
    end
    
    events
        ChangedPoint
    end
    
    % Constructor and destructor
    methods
        function F = filterAndPoint(varargin)
            % set filter
            F.F = xplr.filter(varargin{:});
            if ~ismember(F.F.type,{'selection' 'indices'})
                error 'filterAndPoint filter can be only of type ''selection'' or ''indices'''
            end
            F.headerin = F.F.headerin;
            F.listeners.filter = addlistener(F.F,'ChangedOperation',@(u,e)transitNotification(F,'filter',e));
            
            % set P
            for i=1:F.ndin
                F.P(i) = xplr.point(F.headerin(i));
            end
            for i=1:F.ndin
                F.listeners.point(i) = addlistener(F.P(i),'ChangedPoint',@(u,e)transitNotification(F,'point',e));
                F.listeners.pointselection(i) = addlistener(F.P(i),'ChangedOperation',@(u,e)transitNotification(F,'pointselection',e));
            end
            F.dolistenpoint = true;
            
            % set output header (uses filters or P depending on whether
            % filters selection is non-empty)
            setHeaderout(F)
        end
        function delete(F)
            if ~isvalid(F) && ~isprop(F,'listeners'), return, end
            deleteValid(F.listeners)
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
                    % no tracking
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
        function [headvalue affectedcolumns] = setAddHeaderInfo(F,~,~) %#ok<INUSD,STOUT>
            error 'method ''setAddHeaderInfo'' is not valid for filterAndPoint object; call it rather on the children filter object'
        end
    end
    
    % Point modification
    methods
        function dataind = get.indices(F)
            if F.F.nsel==0
                dataind = F.F.indices;
            else
                dataind = {fn_indices(F.szin,[F.P.index])};
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
            notify(F,'ChangedPoint')
            if ~isequal(F.index,curind) && F.F.nsel==0
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
                        notify(F,'ChangedOperation',xplr.eventinfo('filter','all'))
                    elseif strcmp(e.flag,'new') && ismember(1,e.ind)
                        % new selection(s) replaces point selection by filter
                        % selection
                        notify(F,'ChangedOperation',xplr.eventinfo('filter','all'))
                    else
                        notify(F,'ChangedOperation',e)
                    end
                case 'point'
                    if F.dolistenpoint % it is faster to 'disable' the listener by using this variable than by actually disabling it!!
                        notify(F,'ChangedPoint')
                    end
                case 'pointselection'
                    if F.F.nsel==0 && F.dolistenpoint
                        notify(F,'ChangedOperation',xplr.eventinfo('filter','point'))
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
        function slice = operation(F,x,dims)
            if F.F.nsel>0
                slice = operation(F.F,x,dims);
            else
                % need to do by hand to set the accurate header
                % check input
                checkdata(F,x,dims)
                % slice
                slic = slicing(F.P,x.data,dims,F.ndout);
                % header
                head = x.header;
                ndx = length(head);
                head = [head(1:dims(1)-1) F.headerout head(setdiff(dims(1):ndx,dims))];
                % output
                slice = xplr.xdata(slic,head);
            end
        end
        function updateOperation(F,x,dims,slice,flag,ind)
            % check input
            checkdata(F,x,dims)
            
            % slice
            switch flag
                case 'all'
                    if F.F.nsel>0
                        slic = slicing(F.F,x.data,dims);
                    else
                        if F.ndout~=1, error 'output header should be scalar', end
                        slic = slicing(F.P,x.data,dims,F.ndout);
                    end
                case {'new' 'chg'}
                    slic = slicing(F.F,x.data,dims,ind);
                case 'chg&new'
                    slic = slicing(F.F,x.data,dims,[ind{:}]);
                case 'chg&rm'
                    slic = slicing(F.F,x.data,dims,ind{1});
                case 'remove'
                    if F.F.nsel>0
                        slic = [];
                    else
                        flag = 'all';
                        if F.ndout~=1, error 'output header should be scalar', end
                        slic = slicing(F.P,x.data,dims,F.ndout);
                    end
                case 'perm'
                    slic = [];
                case 'point'
                    if F.ndout~=1, error 'output header should be scalar', end
                    flag = 'chg'; ind = 1;
                    slic = slicing(F.P,x.data,dims,F.ndout);
                otherwise
                    error('flag ''%s'' not handled',flag)
            end
            slice.updateData(flag,dims,ind,slic,F.headerout); % this will trigger automatic notifications
        end
    end
    
end