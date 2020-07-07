classdef DimensionLabel
% DimensionLabel defines a dimension
% function L = DimensionLabel(label, type[, unit|all_units])
% 
% Input:
% * label   a string (e.g. 'time')
% * type    'numeric', 'logical', 'char' or 'mixed'
% * unit    string, cell array with 2 columns, or struct with fields unit and value
%
% Examples:
% tlabel = DimensionLabel('time','numeric',{'s' 1; 'ms' 1e-3; 'min' 60; 'hour' 3600};
% clabel = DimensionLabel('condition','char');   
    
    properties (SetAccess='private')
        label
        type
        unit
        all_units
    end
    properties (Dependent, SetAccess='private', Transient)
        default_val
    end
    
    methods
        function L = DimensionLabel(label, type, unit)
            if ~ischar(label), error 'label must be a character array', end
            L.label = label;
            if ~ismember(type, {'numeric', 'logical', 'char', 'mixed'}), error 'type must be either ''numeric'', ''logical'' or ''char''', end
            L.type = type;
            if nargin<3
                return
            elseif ~strcmp(type, 'numeric')
                error 'unit can be defined only for a ''numeric'' label'
            end
            if ischar(unit)
                L.unit = unit;
                L.all_units = struct('unit', unit, 'value', 1);
            else
                if iscell(unit)
                    unit = cell2struct(unit, {'unit', 'value'}); 
                elseif ~isstruct(unit)
                    error 'invalid definition of unit(s)'
                end
                idxref = find([unit.value] == 1, 1, 'first');
                if isempty(idxref), error 'at least one unit must have value equal to 1', end
                L.unit = unit(idxref).unit;
                [~, ord] = sort([unit.value]);
                L.all_units = unit(ord);
            end
        end
        function x = get.default_val(L)
            x = get_default_value(L.type);
        end
    end
    
    methods (Static)
        function [type, default_val] = infer_type(x)
            if isnumeric(x)
                type = 'numeric';
            elseif islogical(x)
                type = 'logical';
            elseif ischar(x)
                type = 'char';
            else
                type = 'mixed';
            end
            if nargout >= 2, default_val = get_default_value(type); end
        end
    end
    
end


%---
function default_val = get_default_value(type)
    default_val = brick.switch_case(type, ...
        'numeric',  0, ...
        'logical',  false, ...
        'char',     '', ...
        'mixed',    []);
end
