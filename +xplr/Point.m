classdef Point < xplr.DataOperand
   
    properties (AbortSet=true)
        index_exact = 1; % real value
    end
    properties (SetAccess='private')
        index = 1;  % integer between 1 and header_in.n
    end
    properties (Dependent, Transient)
        value  % real-world position
        value_str    % real-world position, with unit
    end
   
    % Constructor and setting index
    methods
        function P = Point(header_in)
            if ~isscalar(header_in)
                P = xplr.Point.empty(1, 0);
                for i=1:length(header_in)
                    P(i) = xplr.Point(header_in(i));
                end
                return
            end
            % no output header because data is averaged to a single value!
            P.header_in = header_in;
            P.header_out = xplr.Header.empty(1, 0);
        end
        function set.index_exact(P, x)
            P.index_exact = x;
            i = max(1, min(P.header_in.n,round(x)));
            chg_ij = (i ~= P.index);
            P.index = i;
            % notification
            notify(P, 'ChangedOperation', xplr.EventInfo('point', chg_ij))
        end
        function copy_in(P, obj)
            P.index = obj.index;
        end
    end
    
    % Conversion to real-world
    methods
        function x = get.value(P)
            head = P.header_in;
            switch head.type
                case 'measure'
                    x = head.start + head.scale*P.index_exact;
                case 'categorical'
                    x = head.values(P.index, :); % cell array
                    if isscalar(x), x = x{1}; end
            end
        end
        function set.value(P,x)
            if ischar(x), P.value_str = x; return, end
            head = P.header_in;
            if head.is_measure
                P.index = (x - head.start)/head.scale;
            else
                P.index = x;
            end
        end
        function str = get.value_str(P)
            x = P.value;
            head = P.header_in;
            switch head.type
                case 'measure'
                    % look for the most appropriate unit!
                    [u, iu] = unique([head.all_units.value]); % if several units have the same value, consider only the first one
                    idx = find(abs(x) >= u, 1, 'last');              % units are oredered by increasing value
                    if isempty(idx), idx = 1; end
                    idx = iu(idx);
                    str = [num2str(x/head.all_units(idx).value), head.all_units(idx).unit];
                case 'categorical'
                    str = x;
            end
        end
        function set.value_str(P, str)
            head = P.header_in;
            switch head.type
                case 'measure'
                    tokens = regexp(str, '^([-\d\.]*)(.*)$', 'tokens');
                    if isempty(tokens), error 'could not read string', end
                    tokens = tokens{1};
                    x = str2double(tokens{1});
                    unit = tokens{2};
                    if ~strcmpi(unit, head.unit)
                        idx = find(strcmpi(unit, {head.all_units.unit}));
                        if isempty(idx), error 'unit is not recognized', end
                        x = x * head.all_units(idx).value;
                    end
                    P.value = x;
                case 'categorical'
                    idx = brick.find(str, head.values, 'rows');
                    if isempty(idx), error 'not a possible value', end
                    P.index = idx;
            end
        end
    end
    
    % Slicing
    methods
        function slic = slicing(P, dat, dims, nd_out)
            % here P can be non-scalar!
            if length(dims) ~= length(P), error 'number of dimensions does not match number of points', end
            if nargin < 4, nd_out = 0; end
            
            % size
            s = size(dat);
            nddata = max(max(dims), length(s));
            s(end + 1:nddata) = 1;
            
            % slice
            subs = substruct('()', repmat({':'}, 1, length(s)));
            for i=1:length(P), subs.subs{dims(i)} = P(i).index; end
            slic = subsref(dat,subs);
            rsh = s;
            switch nd_out % does slicing output space span zero [default] or one dimension?
                case 0
                    rsh(dims) = [];
                case 1
                    rsh(dims(1)) = 1;
                    rsh(dims(2:end)) = [];
                otherwise
                    error 'point slicing output can only occupy zero or one dimension'
            end
            slic = reshape(slic, [rsh, 1]);
        end
    end
    methods (Access='protected')
        function slic = operation_(P, dat, dims)
            % function slic = operation_(P,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            % here again P can be non-scalar...
            slic = slicing(P, dat, dims);
        end
        function update_operation_(P, x, dims, slice)
            slic = slicing(P, x.data, dims);
            slice.chg_data(slic); % this will trigger automatic notifications
        end
    end
    
    % Link with point selection in real world coordinates
    methods
        function point_world = operation_data_to_space(P)
            point_world = P.header_in.start + (P.index_exact-1)*P.header_in.scale;
        end
        function update_operation_data_to_space(P, WO, ~)
            WO.operation = P.operation_data_to_space();
            notify(WO, 'ChangedOperation')
        end
        function update_operation_space_to_data(P, point_world, ~)
            P.index_exact = 1 + (point_world - P.header_in.start)/P.header_in.scale;
        end
    end
    
end
