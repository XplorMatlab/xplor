classdef XData < xplr.GraphNode
    % function x = xdata(dat[,head[,name]])
    %
    % A container for data, associated with header information.
    % A number of different events are thrown when the data or header
    % information is being changed
    % When headers are not provided, opens a graphic interface allowing
    % user to set the headers.
    %
    % Input:
    % * dat   ND array
    % * head  a cell array with as many elements as data dimensions, each element is itself a cell array containing arguments for the xplr.header constructor
    % * name  string
    %
    % See also xplr.header
    
    properties (SetAccess='private')
        data        % ND array
        header = xplr.DimHeader.empty(1, 0);
        name = '';
    end
    
    properties (Dependent, SetAccess='private')
        sz
        nd
    end
    
    events
        ChangedData % sent with info xplr.eventinfo('data',chg_head)
    end
    
    % Constructor and simple access
    methods
        function x = XData(dat, head, name, dim_ids)
            if nargin == 0, return, end % default empty data
            
            % create xplr.header objects
            if nargin < 2 || isempty(head)
                % open user edition window to edit headers
                head = xplr.EditHeader(dat);
            elseif ischar(head)
                head = {head};
            end
            if iscell(head)
                % create xplr.header objects from labels
                labels = head;
                head = xplr.DimHeader.empty(1, 0);
                for i=1:length(labels)
                    if iscell(labels{i})
                        head(i) = xplr.Header(labels{i}{:});
                    else
                        head(i) = xplr.Header(labels{i}, size(dat,i));
                    end
                end
            end
            
            % convert to xplr.DimHeader objects
            if isa(head,'xplr.DimHeader')
                if nargin>=4
                    error 'cannot set dim_ids when header are already dimheader objects'
                end
            else
                if nargin < 4
                    dim_ids = rand(1, length(head));
                end
                head = xplr.DimHeader(head, dim_ids);
            end
            
            % set object
            x.update_data_dim('global', [], dat, head)
            if nargin >= 3
                x.name = name;
            end
        end
        function y = copy(x)
            y = xplr.XData(x.data, x.header);
        end
        function nd = get.nd(x)
            nd = length(x.header);
        end
        function s = get.sz(x)
            s = [x.header.n];
        end
        function [dim, dim_id] = dimension_number_and_id(x, d)
            % function [dim, dim_id] = dimension_number_and_id(x,d)
            %---
            % Convert any of dimension numbers, identifiers or labels
            % to both dimension numbers and identifiers.
            % Returns [] if some dimension was not found.
            [dim, dim_id] = dimension_number_and_id(x.header, d);
        end
        function dim_id = dimension_id(x, d)
            [~, dim_id] = x.dimension_number_and_id(d);
        end
        function dim = dimension_number(x, d)
            [dim, ~] = x.dimension_number_and_id(d);
        end
        function label = dimension_label(x, d)
            label = x.header.dimension_label(d);
        end
        function head = header_by_id(x, dim_id)
            d = x.dimension_number(dim_id);
            head = x.header(d);
        end
        function head = non_singleton_header(x)
            head = x.header(x.sz>1);
        end
    end
    
    % Modification of an xdata object and raising of the corresponding
    % notification
    methods
        function set_name(x, name)
            if strcmp(name, x.name), return, end
            x.name = name;
            notify(x, 'ChangedData', xplr.EventInfo('data', 'name'))
        end
        function chg_Data(x, data)
            if isequal(data, x.data), return, end
            % changes in size are allowed only in the 'measure' dimensions
            data_sz = strictsize(data, x.nd);
            if length(data_sz) > x.nd, error 'Cannot increase number of dimensions. Use update_data to change both data and headers.', end
            chg_sz = (data_sz ~= x.sz);
            if any([x.header(chg_sz).is_categorical_with_values]), error 'Cannot change data size in categorical dimensions. Use update_data to change both data and headers', end
            % set data
            if ~isreal(data), error 'data cannot be complex', end
            x.data = data;
            if ~any(chg_sz)
                notify(x, 'ChangedData', xplr.EventInfo('data', 'chg_data'))
                return
            end
            % update dimension headers if necessary
            for i=find(chg_sz)
                x.header(i) = update_measure_header(x.header(i), data_sz(i));
            end
            notify(x, 'ChangedData', xplr.EventInfo('data', 'chg_dim', find(chg_sz))) %#ok<FNDSB>
        end
        function update_data(x, flag, d, ind, value, new_head)
            % function update_data(x,flag,dim,ind,value,new_head)
            %---
            % arguments value and new_head are supposed to be only the
            % updated parts, i.e.:
            % - for flags 'new' and 'chg', value is the data only for
            %   indices ind
            % - new_head is only the header in dimension dim
            % however giving the full updated data for value, or the full
            % header for new_head, is tolerated
            
            dim = x.dimension_number(d);
            if isempty(dim), error('dimension %g absent from data', d); end
            
            % check that value is real
            if nargin >= 5 && ~isreal(value), error 'data cannot be complex', end
            
            % check that header is a dimheader
            if ~isa(new_head, 'xplr.DimHeader'), error 'new header must be a dimheader object', end
            
            % update header
            if nargin>=6
                if length(new_head) == x.nd
                    % giving the full headers instead of only the updated one
                    new_head = new_head(dim);
                end
                check_header_update(x.header(dim), flag, ind, new_head)
            else
                new_head = update_header(x.header(dim), flag, ind);
            end
            if strcmp(flag, 'all') && isequal(new_head, x.header(dim))
                % flag 'chg_data' might be preferable to 'all' to indicate
                % that header did not change
                flag = 'chg_data';
            end
            tmp = new_head;
            new_head = x.header;
            new_head(dim) = tmp;
            clear tmp
            % update data
            if ~fn_ismemberstr(flag, {'chg_data', 'all', 'chg_dim'})
                s = substruct('()', repmat({':'}, 1, x.nd));
                s.subs{dim} = ind;
            end
            switch flag
                case {'all', 'chg_data', 'chg_dim'}
                    new_data = value;
                case {'chg', 'new', 'chg&new', 'chg&rm'}
                    if size(value,dim) == new_head(dim).n
                        % giving the full data instead of only the updated part
                        new_data = value;
                    else
                        if strcmp(flag, 'chg&new')
                            s.subs{dim} = [ind{:}];
                        elseif strcmp(flag, 'chg&rm')
                            s.subs{dim} = ind{1};
                        end
                        new_data = subsasgn(x.data, s, value);
                        if strcmp(flag, 'chg&rm')
                            s.subs{dim} = ind{2};
                            new_data = subsasgn(new_data, s, []);
                        end
                    end
                case 'remove'
                    new_data = subsasgn(x.data, s, []);
                case 'perm'
                    new_data = subsref(x.data, s);
                otherwise
                    error('invalid flag ''%s'' for xdata update_data method')
            end
            if ~isequal(strictsize(new_data, length(new_head)), [new_head.n])
                error 'new data size does not match header(s)'
            end
            % really update data only now (after all checks occured)
            x.header = new_head;
            x.data = new_data;
            % notification
            notify(x, 'ChangedData', xplr.EventInfo('data', flag, dim, ind))
        end
        function update_data_dim(x, flag, dim, new_data, new_head)
            % update header
            if ~isa(new_head,'xplr.DimHeader'), error 'new header must be a dimheader object', end
            switch flag
                case 'global'
                    x.header = new_head;
                case 'chg_dim'
                    if length(dim) ~= length(new_head), error 'length of new header does not match number of new dimensions', end
                    x.header(dim) = new_head;
                otherwise
                    error('invalid flag ''%s'' for xdata update_data_dim method')
            end
            % update data
            if ~isreal(new_data), error 'data cannot be complex', end
            if x.nd == 1 && isvector(new_data), new_data = new_data(:); end % don't generate an error for a row vector input
            if ~isequal(strictsize(new_data, x.nd), x.sz)
                error 'new data size does not match new header'
            end
            x.data = new_data;
            % notification
            notify(x, 'ChangedData', xplr.EventInfo('data', flag, [x.header(dim).dim_id]))
        end
    end
    
end
