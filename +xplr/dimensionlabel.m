classdef dimensionlabel
    % function L = dimensionlabel(label,type[,unit|allunits])
    %---
    % Define a dimension.
    % 
    % Input:
    % - label   a string (e.g. 'time')
    % - type    'numeric', 'logical', 'char' or 'object'
    % - unit    string, cell array with 2 columns, or struct with fields
    %           unit and value
    %
    % Examples:
    % tlabel = dimensionlabel('time','numeric',{'s' 1; 'ms' 1e-3; 'min' 60; 'hour' 3600};
    % clabel = dimensionlabel('condition','char');   
    
    properties (SetAccess='private')
        label
        type
        unit
        allunits
    end
    properties (Dependent, SetAccess='private')
        defaultval
    end
    
    methods
        function L = dimensionlabel(label,type,unit)
            if ~ischar(label), error 'label must be a character array', end
            L.label = label;
            if ~ismember(type,{'numeric' 'logical' 'char' 'object'}), error 'type must be either ''numeric'', ''logical'' or ''char''', end
            L.type = type;
            if nargin<3
                return
            elseif ~strcmp(type,'numeric')
                error 'unit can be defined only for a ''numeric'' label'
            end
            if ischar(unit)
                L.unit = unit;
                L.allunits = struct('unit',unit,'value',1);
            else
                if iscell(unit)
                    unit = cell2struct(unit,{'unit' 'value'}); 
                elseif ~isstruct(unit)
                    error 'invalid definition of unit(s)'
                end
                idxref = find([unit.value]==1,1,'first');
                if isempty(idxref), error 'at least one unit must have value equal to 1', end
                L.unit = unit(idxref).unit;
                [~, ord] = sort([unit.value]);
                L.allunits = unit(ord);
            end
        end
        function x = get.defaultval(L)
            x = getDefaultValue(L.type);
        end
    end
    
    methods (Static)
        function [type defaultval] = infertype(x)
            if isnumeric(x)
                type = 'numeric';
            elseif islogical(x)
                type = 'logical';
            elseif ischar(x)
                type = 'char';
            else
                error('header values cannot be of class ''%s''',class(x))
            end
            if nargout>=2, defaultval = getDefaultValue(type); end
        end
    end
    
end


%---
function defaultval = getDefaultValue(type)
    if strcmp(type,'object'), error 'no default value for ''object'' type', end
    defaultval = fn_switch(type, ...
        'numeric',  0, ...
        'logical',  false, ...
        'char',     '');
end

