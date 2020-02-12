classdef eventinfo < event.EventData & dynamicprops
    % function obj = eventinfo(type,arguments...)
    % type: 'filter'
    % arguments:
    % * 'all'
    % * 'new',ind
    % * 'chg',ind
    % * 'remove',ind
    % * 'chg&new',ind
    % * 'chg&rm',ind
    % * 'perm',ind
    % * 'point'
    %
    % type: 'point'
    % arguments:    chgij [default false]
    %
    % type: 'data'
    % arguments:
    % * 'global'
    % * 'chgdata'
    % * 'insertdim',dim
    % * 'rmdim',dim
    % * 'chgdim',dim
    % * 'all',dim
    % * 'new',dim,ind
    % * 'chg',dim,ind
    % * 'remove',dim,ind
    % * 'chg&new',{indchg indnew}
    % * 'chg&rm',{indchg indrm}
    % * 'perm',dim,ind
    %
    % type: 'zoom'
    % arguments:    chgnout,dim
    % 
    % type: 'clip'
    % arguments:
    % * 'clip',value
    % * 'automode'
    % * 'adjust'
    % * 'span'

    properties 
        type
    end
    methods
        function obj = eventinfo(type,varargin)
            obj.type = type;
            switch type
                case 'point'
                    F = {'chgij'};
                case 'data'
                    F = {'flag' 'dim' 'ind'};
                case 'filter'
                    F = {'flag' 'ind' 'value'};
                case 'zoom'
                    F = {'chgnout' 'dim'};
                case 'clip'
                    F = {'flag' 'value'};
                otherwise
                    error('unknown event type ''%s''',type)
            end
            for i=1:length(varargin)
                addprop(obj,F{i});
                obj.(F{i})=varargin{i};
            end
        end
    end
end