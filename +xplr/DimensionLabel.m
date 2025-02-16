classdef DimensionLabel
% DimensionLabel defines a dimension
% function L = DimensionLabel(label, type[, unit|all_units])
% 
% Input:
% * label   a string (e.g. 'time')
% * type    'numeric', 'logical', 'char', 'color' or 'mixed'
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
            assert(ismember(type, {'numeric', 'logical', 'char', 'color', ...
                'datetime', 'duration', ...
                'Selection1D', 'Selection2D', 'mixed'}))
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
            elseif ismember(class(x), {'logical', 'char', 'datetime', 'duration'})
                type = class(x);
            elseif isa(x, 'xplr.SelectionND')
                switch x.nd
                    case 1
                        type = 'Selection1D';
                    case 2
                        type = 'Selection2D';
                    otherwise
                        error 'SelectionND header values with nd>2 not handled'
                end
            else
                type = 'mixed';
            end
            if nargout >= 2, default_val = get_default_value(type); end
        end
    end
    
end


%---
function default_val = get_default_value(type)
    switch type
        case 'numeric'
            default_val = 0;
        case 'logical'
            default_val = false;
        case 'char'
            default_val = '';
        case 'color'
            default_val = [0, 0, 0];
        case 'Selection1D'
            default_val = xplr.SelectionND('empty', 1);
        case 'Selection2D'
            default_val = xplr.SelectionND('empty', 2);
        case 'mixed'
            default_val = [];
    end
end
