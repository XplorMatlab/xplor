classdef filterAndPoint < xplr.dataoperand
    % function F = filter(headerin[,label])
    %---
    % The class filerAndPoint combines a filter (operating on a
    % n-dimensional space) and a set of n points. 
    
    properties
        F
        P = xplr.point.empty(1,0);
    end
    properties (SetAccess='private')
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
            F.F = xplr.filter(varargin{:});
            F.headerin = F.F.headerin;
            connectlistener(F.F,F,'ChangedOperation',@(u,e)transitNotification(F,'filter',e))
            
            % set P
            for i=1:F.ndin
                F.P(i) = xplr.point(F.headerin(i));
            end
            for i=1:F.ndin
                connectlistener(F.P(i),F,'ChangedPoint',@(u,e)transitNotification(F,'point',e))
                connectlistener(F.P(i),F,'ChangedOperation',@(u,e)transitNotification(F,'pointselection',e))
            end
            F.dolistenpoint = true;
            
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
                    % no tracking: we keep the current value for
                    % F.headerout
                    return
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
                head(dims) = [];
                head = [head(1:dims(1)-1) F.headerout head(dims(1):end)];
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