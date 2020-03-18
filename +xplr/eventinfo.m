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
    % * 'global'                    data is potentially completely different
    % * 'chgdata'                   data has changed but sizes and all header information remain the same
    % * 'chgdim',dim                data has changed, header and size in dimension dim have changed (dim can be non-scalar here, but not in the other options below)
    % * 'all',dim                   data has changed, size and header value tables in dimension dim have changed, but not the header name (for example 'time')
    % * 'new',dim,ind               new data (and therefore header values) have been inserted along dimension dim
    % * 'chg',dim,ind               data has changed (and therefore also header values) at some specific positions in dimension dim
    % * 'remove',dim,ind            some data (and corresdponding header values) have been removed along dimension dim
    % * 'chg&new',{indchg indnew}   some data/header values have been removed and some added in dimension dim (less removals than additions -> seen as change + addition)
    % * 'chg&rm',{indchg indrm}     some data/header values have been removed and some added in dimension dim (more removals than additions -> seen as change + removals)
    % * 'perm',dim,ind              data has not changed but data and header values were permuted along dimension dim
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
            for i=1:length(F)
                addprop(obj,F{i});
            end
            for i=1:length(varargin)
                obj.(F{i})=varargin{i};
            end
        end
    end
end