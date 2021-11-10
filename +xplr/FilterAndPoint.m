classdef FilterAndPoint < xplr.DataOperand
    % function F = filter(header_in[,label])
    % function F = filter(F,P)
    %---
    % The class filerAndPoint combines a filter (operating on a
    % n-dimensional space) and a set of n points. 
    
    properties
        F
        P = xplr.Point.empty(1, 0);
    end
    properties (SetAccess='protected', Transient)
        do_listen_point = true;
    end
    properties (Dependent, SetAccess='protected', AbortSet=true, Transient)
        indices             % current selection
        point_index_exact   % current point
    end
    properties (Dependent, Transient, SetAccess='private')
        point_index         % current point
    end

    % Constructor and destructor
    methods
        function F = FilterAndPoint(varargin)
            % init array of filters?
            if nargin == 0
                return
            elseif isscalar(varargin) && isnumeric(varargin{1})
                n = varargin{1};
                if n == 0
                    F = xplr.FilterAndPoint.empty(1,0);
                else
                    F(n) = xplr.FilterAndPoint();
                end
                return
            end
            
            % set filter
            create_members = ~isa(varargin{1}, 'xplr.Filter');
            if create_members
                F.F = xplr.Filter(varargin{:});
            else
                F.F = varargin{1};
            end
            F.header_in = F.F.header_in;
            brick.connect_listener(F.F, F, 'ChangedOperation', @(u,e)transit_notification(F, 'filter', e))
            
            % set P
            if create_members
                for i=1:F.nd_in
                    F.P(i) = xplr.Point(F.header_in(i));
                end
            else
                F.P = varargin{2};
                if length(F.P) ~= F.nd_in, error 'number of point filters does not match number of input dimensions', end
                for i=1:F.nd_in
                    if ~isequal(F.P(i).header_in, F.header_in(i))
                        error 'headers of filter and point filter(s) do not match'
                    end
                end
            end
            for i=1:F.nd_in
                brick.connect_listener(F.P(i), F, 'ChangedOperation', @(u,e)transit_notification(F, 'point', e))
            end
            
            % set output header (uses filters or P depend_ing on whether
            % filters selection is non-empty)
            set_header_out(F)
        end
    end
    
    % Header
    methods
        function set_header_out(F)
            F.header_out = F.F.header_out;
            if F.F.n_sel == 0
                % set header values
                % (tracking of input header values)
                if ~isscalar(F.header_in) 
                    % do not track values in tables of each input dimension
                    % if any
                    n_colin = 1;
                    head_values = {brick.strcat({F.P.index}, ',')};
                elseif F.header_in.n_column > 0
                    head_values = F.header_in.values(F.P.index, :);
                    n_colin = size(head_values, 2);
                elseif F.header_in.categorical
                    n_colin = 1;
                    head_values = {F.P.index};
                else
                    n_colin = 1;
                    head_values = F.header_in.get_item_names(F.P.index);
                end
                % (point selection cannot define additional values: put
                % default values if such additional values were defined
                % during earlier ROI selection)
                for i = n_colin + 1:F.header_out.n_column
                    head_values{1, i} = F.header_out.sub_labels(i).default_val;
                end
                % replace output header
                F.header_out = update_header(F.header_out, 'new', 1, head_values);
            end
        end
        function augment_header(F,varargin)
            augment_header(F.F, varargin{:})
            set_header_out(F)
        end
        function [head_value, affected_columns] = set_add_header_info(F,~,~) %#ok<INUSD,STOUT>
            error 'method ''set_add_header_info'' is not valid for FilterAndPoint object; call it rather on the children filter object'
        end
    end
    
    % Point modification
    methods
        function data_ind = get.indices(F)
            if F.F.n_sel == 0
                data_ind = {brick.indices(F.szin, [F.P.index])};
            else
                data_ind = F.F.indices;
            end
        end
        function x = get.point_index_exact(F)
            x = [F.P.index_exact];
        end
        function x = get.point_index(F)
            x = [F.P.index];
        end
        function set.point_index_exact(F, x)
            
            % if there are multiple points (i.e. input space is
            % multi-dimensional), it is better to generate a single event
            % after all points have been modified rather than transit one
            % event per modified point
            F.do_listen_point = false;
            c = onCleanup(@()set(F, 'do_listen_point', true));
            
            % update point
            cur_ind = F.point_index;
            for i=length(F.P), F.P(i).value = x(i); end
            
            % update header
            set_header_out(F)
            
            % generation of a single event
            chg_ij = ~isequal(F.point_index, cur_ind);
            notify(F, 'ChangedOperation', xplr.EventInfo('point', chg_ij))
            if chg_ij && F.F.n_sel == 0
                notify(F, 'ChangedOperation', xplr.EventInfo('filter', 'chg', 1))
            end
        end
    end
    
    % Selection modification 
    methods
        function update_selection(F, varargin)
            update_selection(F.F, varargin{:})
            set_header_out(F)
        end
    end
    
    % Notification
    methods
        function transit_notification(F, flag, e)
            set_header_out(F)
            switch flag
                case 'filter'
                    if strcmp(e.type, 'operation')
                        if F.F.n_sel == 0
                            return
                        else
                            e2 = xplr.EventInfo('filter', 'chg', 1:F.F.n_sel);
                        end
                    elseif strcmp(e.flag, 'remove') && F.F.n_sel == 0
                        % removal of selection(s) replaces filter selection
                        % by point selection
                        e2 = xplr.EventInfo('filter', 'all');
                    elseif strcmp(e.flag, 'new') && ismember(1, e.ind)
                        % new selection(s) replaces point selection by filter
                        % selection
                        e2 = xplr.EventInfo('filter', 'all');
                    else
                        % avoid very weird Matlab bug when simply passing
                        % e: make a copy of it
                        e2 = xplr.EventInfo('filter', e.flag,e.ind, e.value);
                    end
                    notify(F, 'ChangedOperation', e2)
                case 'point'
                    if F.do_listen_point % 'do_listen_point' property is manipulated in F.set.point_index
                        % same as above: copy e
                        e2 = xplr.EventInfo('point', e.chg_ij);
                        notify(F, 'ChangedOperation', e2)
                        if e.chg_ij && F.F.n_sel == 0
                            notify(F, 'ChangedOperation', xplr.EventInfo('filter', 'chg', 1))
                        end
                    end
            end
        end
    end
    
    % Slicing
    methods
        function slic = slicing(F, dat, dims, sel_sub_idx)
            if F.F.n_sel > 0
                if nargin < 4, sel_sub_idx = 1:F.F.n_sel; end
                slic = slicing(F.F, dat, dims, sel_sub_idx);
            else
                if nargin >= 4 && ~isequal(sel_sub_idx, 1), error 'invalid selection indices', end
                slic = slicing(F.P, dat, dims, F.nd_out);
            end
        end
    end
    methods (Access='protected')
        function slic = operation_(F, dat, dims)
            % function slic = operation_(F,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            if F.F.n_sel > 0
                slic = F.F.operation_(dat, dims);
            else
                nd_out = 1; % specify that slicing dimension should not be removed
                slic = F.P.slicing(dat, dims, nd_out);
            end
        end
        function update_operation_(F, x, dims, slice, flag, ind)
            if F.F.n_sel > 0
                F.F.update_operation_(F, x, dims, slice, flag, ind)
            else
                F.P.update_operation_(F, x, dims, slice)
            end
            %             % slice
            %             switch flag
            %                 case 'all'
            %                     if F.F.n_sel>0
            %                         slic = slicing(F.F,x.data,dims);
            %                     else
            %                         if F.nd_out~=1, error 'output header should be scalar', end
            %                         slic = slicing(F.P,x.data,dims,F.nd_out);
            %                     end
            %                 case {'new' 'chg'}
            %                     slic = slicing(F.F,x.data,dims,ind);
            %                 case 'chg&new'
            %                     slic = slicing(F.F,x.data,dims,[ind{:}]);
            %                 case 'chg&rm'
            %                     slic = slicing(F.F,x.data,dims,ind{1});
            %                 case 'remove'
            %                     if F.F.n_sel>0
            %                         slic = [];
            %                     else
            %                         flag = 'all';
            %                         if F.nd_out~=1, error 'output header should be scalar', end
            %                         slic = slicing(F.P,x.data,dims,F.nd_out);
            %                     end
            %                 case 'perm'
            %                     slic = [];
            %                 otherwise
            %                     error('flag ''%s'' not handled',flag)
            %             end
            %             slice.updateData(flag,dims,ind,slic,F.header_out); % this will trigger automatic notifications
        end
    end
    
    % Link with selection in real world coordinates
    methods
        function selection_world = operation_data_to_space(F)
            error 'FilterAndPoint object should not be directly connected to a wordOperand object: connect rather its child point and filter'
        end
        function update_operation_data_to_space(F, WO, e)
            error 'FilterAndPoint object should not be directly connected to a wordOperand object: connect rather its child point and filter'
        end
        function update_operation_space_to_data(F, ~, e)
            error 'FilterAndPoint object should not be directly connected to a wordOperand object: connect rather its child point and filter'
        end
    end
   
    % Copy (see xplr.dataOperand.loadfromfile)
    methods
        function copy_in(F, obj)
            F.F.copy_in(obj.F)
            F.P.copy_in(obj.P)
        end
    end
    
    % Context menu
    methods
        function context_menu(F, m)
            % the context menu of the FilterAndPoint object is no more than
            % the context menu of its filter component object
            context_menu(F.F, m)
        end
    end
end
