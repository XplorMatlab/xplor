classdef bankRegistry < handle
% bank registry
%
    properties (Dependent, SetAccess='private')
        n
    end
    properties
        content = struct('key',cell(1,0),'value',cell(1,0),'users',cell(1,0));
    end
    
    methods
        function n = get.n(R)
            n = length(R.content);
        end
    end
    
    methods (Access='private')
        function idx = getIndex(R,key)
            idx = fn_find(key,{R.content.key},'first');
        end
    end
    methods
        function register(R,key,value,user)
            douser = (nargin>=4);
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = getValue(R,key{1});
                    if isempty(R1)
                        % create sub-register
                        R1 = xplr.bankRegistry();
                        idx = R.n+1;
                        R.content(idx) = struct('key',key{1},'value',R1,'users',{cell(1,0)});
                    elseif ~isa(R1,'xplr.bankRegistry')
                        error 'key is invalid: encountered a leave instead of a sub-registry'
                    end
                    if douser
                        register(R1,key(2:end),value,user)
                    else
                        register(R1,key(2:end),value)
                    end
                    return
                end
            end
            % Verify that no entry exists yet
            if ~isempty(getIndex(R,key))
                if ischar(key) || isnumeric(key)
                    error('Register already has a value for key ''%s''.',num2str(key))
                else
                    error('Register already has a value for the specified key.')
                end
            end
            % Create new entry
            idx = R.n+1;
            if douser
                R.content(idx) = struct('key',key,'value',value,'users',{{user}});
                if isobject(user)
                    addlistener(user,'ObjectBeingDestroyed', ...
                        @(u,e)R.unregister(key,user));
                end
                xplr.debuginfo('registry', 'key %s 1st user %s (value %s)', ...
                    fn_hash(key,3), user, value)
            else
                R.content(idx) = struct('key',key,'value',value,'users',{cell(1,0)});
                xplr.debuginfo('registry', 'key %s no user (value %s)', ...
                	fn_hash(key,3), value)
            end
        end
        function removed = unregister(R,key,user)
            if ~isvalid(R), return, end % R has been deleted already
            douser = (nargin>=3);
            rm = [];
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = getValue(R,key{1});
                    if isempty(R1)
                        disp 'sub-registry where to unregister was not found'
                        if nargout, removed = false; end
                        return
                    elseif ~isa(R1,'xplr.bankRegistry')
                        error 'key is invalid: encoutered a leave instead of a sub-registry'
                    end
                    if douser
                        rm = unregister(R1,key(2:end),user);
                    else
                        rm = unregister(R1,key(2:end));
                    end
                    if R1.n==0
                        % sub-registry still has no more entries, will be
                        % removed below
                        key = key{1};
                        douser = false;
                    else
                        % sub-registry still has entries, do not remove!
                        if nargout, removed = rm; end
                        return
                    end
                end
            end
            % Unregister if no more users
            idx = getIndex(R,key);
            if isempty(idx)
                disp 'entry to unregister was not found'
                if isempty(rm), rm = false; end
                if nargout, removed = rm; end
                return
            end
            if douser
                idxuser = fn_find(user,R.content(idx).users,'first');
                if isempty(idxuser)
                    disp 'user to unregister was not found in users list'
                else
                    c = R.content(idx);
                    xplr.debuginfo('registry', 'key %s rmv user %s', ...
                        fn_hash(c.key,3), c.users{idxuser})
                    R.content(idx).users(idxuser) = [];                    
                end
            end
            if isempty(R.content(idx).users)
                xplr.debuginfo('registry', 'key %s remove key', ...
                    fn_hash(R.content(idx).key,3))
                R.content(idx) = [];
                if isempty(rm), rm = true; end
            else
                if isempty(rm), rm = false; end
            end
            % Output?
            if nargout, removed = rm; end
        end
        function value = getValue(R,key,newuser)
            douser = (nargin>=3);
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = getValue(R,key{1});
                    if isempty(R1)
                        disp 'sub-registry was not found'
                        return
                    elseif ~isa(R1,'xplr.bankRegistry')
                        error 'key is invalid: encoutered a leave instead of a sub-registry'
                    end
                    if douser
                        value = getValue(R1,key(2:end),newuser);
                    else
                        value = getValue(R1,key(2:end));
                    end
                    return
                end
            end
            % Get value
            idx = getIndex(R,key);
            if isempty(idx)
                value = [];
                return
            end
            value = R.content(idx).value;
            % Add new user
            if douser
                idxuser = fn_find(newuser,R.content(idx).users,'first');
                if isempty(idxuser)
                    xplr.debuginfo('registry', 'key %s add user %s', ...
                        fn_hash(R.content(idx).key,3), newuser)
                    R.content(idx).users{end+1} = newuser;
                    if isobject(newuser)
                        addlistener(newuser,'ObjectBeingDestroyed', ...
                            @(u,e)R.unregister(key,newuser));
                    end
                else
                    disp 'new user to register was already found in users list'
                end
            end
        end
        function clear(R)
            R.content(:) = [];
        end
    end
    
end