classdef Parameters < handle
    % XPLR.PARAMETERS  handle parameters stored in a xml file
    % has the following methods
    % xplr.parameters.get_all_par()  get all parameters
    % xplr.parameters.reload()     reload from file
    % xplr.parameters.get(key)     get a specific parameter
    % xplr.parameters.set(key, value)  set value of a parameter
    
    properties
        params
    end
   
    % Constructor is private
    methods (Access='private')
        function P = Parameters
        end
    end
    
    % Only static functions are public
    methods (Static)
        function P = get_all_par(force_reload)
            persistent P_mem
            if nargin < 1, force_reload = false; end
            if isempty(P_mem) || force_reload
                f_name = fullfile(fileparts(which('xplor')), 'xplor parameters.xml');
                if exist(f_name, 'file')
                    s = fn_readxml(f_name);
                else
                    s = struct;
                end
                P_mem = xplr.Parameters;
                P_mem.params = s;
            end
            P = P_mem.params;
        end
        function reload()
            xplr.Parameters.get_all_par(true);
        end
        function value = get(str)
            value = xplr.Parameters.get_all_par();
            if nargin
                strc = fn_strcut(str, '.');
                for i=1:length(strc)
                    value = value.(strc{i});
                end
            end
        end
        function set(str, value)
            % check value
            if isnumeric(value) || islogical(value)
                if ~isscalar(value), error 'numerical or logical values must be scalar', end
            elseif ~ischar(value)
                error 'value is not a valid parameter'
            end
            % get parameter structure
            s = xplr.Parameters.get_all_par();
            % set value
            str = fn_strcut(str, '.');
            s = set_struct(s, str, value);
            % save
            f_name = fullfile(fileparts(which('xplor')), 'xplor parameters.xml');
            fn_savexml(f_name, s)
            % reload
            xplr.Parameters.get_all_par(true);
        end
    end
    
end


%---
function s = set_struct(s, str, value)

    if isscalar(str)
        s.(str{1}) = value;
    else
        if isfield(s, str{1})
            s1 = s.(str{1});
        else
            s1 = struct;
        end
        s.(str{1}) = set_struct(s1, str(2:end), value);
    end

end
