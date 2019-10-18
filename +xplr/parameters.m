classdef parameters < handle
    % XPLR.PARAMETERS  handle parameters stored in a xml file
    % has the following methods
    % xplr.getAllPar()  get all parameters
    % xplr.reload()     reload from file
    % xplr.get(key)     get a specific parameter
    % xplr.set(key, value)  set value of a parameter
    
    properties
        params
    end
   
    % Constructor is private
    methods (Access='private')
        function P= parameters
        end
    end
    
    % Only static functions are public
    methods (Static)
        function P = getAllPar(forcereload)
            persistent Pmem
            if nargin<1, forcereload = false; end
            if isempty(Pmem) || forcereload
                fname = fullfile(fileparts(which('xplor')),'xplor parameters.xml');
                if exist(fname,'file')
                    s = fn_readxml(fname);
                else
                    s = struct;
                end
                Pmem = xplr.parameters;
                Pmem.params = s;
            end
            P = Pmem.params;
        end
        function reload()
            xplr.parameters.getAllPar(true);
        end
        function value = get(str)
            value = xplr.parameters.getAllPar();
            if nargin
                strc = fn_strcut(str,'.');
                for i=1:length(strc)
                    value = value.(strc{i});
                end
            end
        end
        function set(str,value)
            % check value
            if isnumeric(value) || islogical(value)
                if ~isscalar(value), error 'numerical or logical values must be scalar', end
            elseif ~ischar(value)
                error 'value is not a valid parameter'
            end
            % get parameter structure
            s = xplr.parameters.getAllPar();
            % set value
            str = fn_strcut(str,'.');
            s = setstruct(s,str,value);
            % save
            P.params = s;
            fname = fullfile(fileparts(which('xplor')),'xplor parameters.xml');
            fn_savexml(fname,s)
        end
    end
    
end


%---
function s = setstruct(s,str,value)

    if isscalar(str)
        s.(str{1}) = value;
    else
        if isfield(s,str{1})
            s1 = s.(str{1});
        else
            s1 = struct;
        end
        s.(str{1}) = setstruct(s1,str(2:end),value);
    end

end

