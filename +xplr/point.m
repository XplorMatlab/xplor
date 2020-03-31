classdef point < xplr.dataOperand
   
    properties (SetAccess='private')
        index0 = 1; % real value
    end
    properties
        index = 1;  % integer between 1 and headerin.n
    end
    properties (Dependent)
        value  % real-world position
        valuestr    % real-world position, with unit
    end
   
    % Constructor and setting index
    methods
        function P = point(headerin)
            if ~isscalar(headerin)
                P = xplr.point.empty(1,0);
                for i=1:length(headerin)
                    P(i) = xplr.point(headerin(i));
                end
                return
            end
            % no output header because data is averaged to a single value!
            P.headerin = headerin;
            P.headerout = xplr.header.empty(1,0);
        end
        function set.index(P,x)
            if x==P.index0, return, end %#ok<*MCSUP>
            P.index0 = x;
            i = max(1,min(P.headerin.n,round(x)));
            chgij = (i~=P.index);
            P.index = i;
            % notification
            notify(P,'ChangedOperation',xplr.eventinfo('point',chgij))
        end
    end
    
    % Conversion to real-world
    methods
        function x = get.value(P)
            head = P.headerin;
            switch head.type
                case 'measure'
                    x = head.start + head.scale*P.index0;
                case 'categorical'
                    x = head.values(P.index,:); % cell array
                    if isscalar(x), x = x{1}; end
            end
        end
        function set.value(P,x)
            if ischar(x), P.valuestr = x; return, end
            head = P.headerin;
            if head.ismeasure
                P.index = (x-head.start)/head.scale;
            else
                P.index = x;
            end
        end
        function str = get.valuestr(P)
            x = P.value;
            head = P.headerin;
            switch head.type
                case 'measure'
                    % look for the most appropriate unit!
                    [u iu] = unique([head.allunits.value]); % if several units have the same value, consider only the first one
                    idx = find(abs(x)>=u,1,'last');              % units are oredered by increasing value
                    if isempty(idx), idx = 1; end
                    idx = iu(idx);
                    str = [num2str(x/head.allunits(idx).value) head.allunits(idx).unit];
                case 'categorical'
                    str = x;
            end
        end
        function set.valuestr(P,str)
            head = P.headerin;
            switch head.type
                case 'measure'
                    tokens = regexp(str,'^([-\d\.]*)(.*)$','tokens');
                    if isempty(tokens), error 'could not read string', end
                    tokens = tokens{1};
                    x = str2double(tokens{1}); unit = tokens{2};
                    if ~strcmpi(unit,head.unit)
                        idx = find(strcmpi(unit,{head.allunits.unit}));
                        if isempty(idx), error 'unit is not recognized', end
                        x = x * head.allunits(idx).value;
                    end
                    P.value = x;
                case 'categorical'
                    idx = fn_find(str,head.values,'rows');
                    if isempty(idx), error 'not a possible value', end
                    P.index = idx;
            end
        end
    end
    
    % Slicing
    methods
        function slic = slicing(P,dat,dims,ndout)
            % here P can be non-scalar!
            if length(dims)~=length(P), error 'number of dimensions does not match number of points', end
            if nargin<4, ndout = 0; end
            
            % size
            s = size(dat);
            nddata = max(max(dims),length(s));
            s(end+1:nddata) = 1;
            
            % slice
            subs = substruct('()',repmat({':'},1,length(s)));
            for i=1:length(P), subs.subs{dims(i)} = P(i).index; end
            slic = subsref(dat,subs);
            rsh = s;
            switch ndout % does slicing output space span zero [default] or one dimension?
                case 0
                    rsh(dims) = [];
                case 1
                    rsh(dims(1)) = 1;
                    rsh(dims(2:end)) = [];
                otherwise
                    error 'point slicing output can only occupy zero or one dimension'
            end
            slic = reshape(slic,[rsh 1]);
        end
    end
    methods (Access='protected')
        function slic = operation_(P,dat,dims)
            % function slic = operation_(P,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            % here again P can be non-scalar...
            slic = slicing(P,dat,dims);
        end
        function updateOperation_(P,x,dims,slice)
            slic = slicing(P,x.data,dims);
            slice.chgData(slic); % this will trigger automatic notifications
        end
    end
    
    % Link with point selection in real world coordinates
    methods
        function point_world = operationData2Space(P)
            point_world = P.headerin.start + (P.index-1)*P.headerin.scale;
        end
        function updateOperationData2Space(P,WO,~)
            WO.operation = P.operationData2Space();
            notify(WO,'ChangedOperation')
        end
        function updateOperationSpace2Data(P,point_world,~)
            P.index = 1 + (point_world - P.headerin.start)/P.headerin.scale;
        end
    end
    
end