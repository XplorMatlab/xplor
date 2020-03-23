classdef bankRegistry < handle
% Bank registry stores values indexed by keys. If key is a cell array, a
% nested registry structure will be used (the first element of the cell
% array is used as a first key, the value is itself a registry indexed by
% the second element of the cell array, etc.). 
% When an object is registered, a user of this object must be passed as an
% additional argument: this will serve to determine when the object can be
% unregistered and deleted, i.e. as soon as all its users either have 
% initiated an unregistration, or have been deleted.

    properties (Dependent, SetAccess='private')
        n
    end
    properties
        content = struct('key',cell(1,0),'value',cell(1,0),'users',cell(1,0));
    end
    
    events
        Empty
    end
    
    methods
        function n = get.n(R)
            n = length(R.content);
        end
        function str = char(R)
            str = 'bankRegistry';
        end
    end
    
    methods (Access='private')
        function idx = getIndex(R,key)
            idx = fn_find(key,{R.content.key},'first');
        end
    end
    methods
        function register(R,key,value,user)
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
                        % watch for it to become empty
                        addlistener(R1,'Empty',@(u,e)unregister(R,key{1}));
                    elseif ~isa(R1,'xplr.bankRegistry')
                        error 'key is invalid: encountered a leave instead of a sub-registry'
                    end
                    register(R1,key(2:end),value,user)
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
            if isobject(user)
                hld = addlistener(user,'ObjectBeingDestroyed', @(u,e)R.unregister(key,user));
            else
                hld = [];
            end
            R.content(idx) = struct('key',key,'value',value,'users',{{user; hld}});
            xplr.debuginfo('registry', 'key %s 1st user %s (value %s)', ...
                fn_hash(key,3), char(user), char(value))
        end
        function removed = unregister(R,key,user)
            if ~isvalid(R), return, end % R has been deleted already
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
                    if nargin<3
                        rm = unregister(R1,key(2:end)); % note that if R1 becomes empty, it will emit 'Empty' and automatically get removed from R and be deleted
                    else
                        rm = unregister(R1,key(2:end),user); % note that if R1 becomes empty, it will emit 'Empty' and automatically get removed from R and be deleted
                    end
                    if nargout, removed = rm; end
                    return
                end
            end
            % Value index
            idx = getIndex(R,key);
            if isempty(idx)
                xplr.debuginfo('registry', 'key %s was not found', fn_hash(key,3))
                disp 'entry to unregister was not found'
                if nargout, removed = false; end
                return
            end
            % Unregister value if no more users
            if nargin<3
                % remove value regardless of there are users or not
                R.content(idx).users = [];
            else
                idxuser = fn_find(user,R.content(idx).users(1,:),'first');
                if isempty(idxuser)
                    disp 'user to unregister was not found in users list'
                else
                    c = R.content(idx);
                    xplr.debuginfo('registry', 'key %s rmv user %s', ...
                        fn_hash(c.key,3), char(c.users{1,idxuser}))
                    hld = R.content(idx).users{2,idxuser}; % handle of listener listening for object deletion
                    delete(hld)
                    R.content(idx).users(:,idxuser) = [];                    
                end
            end
            if isempty(R.content(idx).users)
                % no more user -> unregister value and delete if (if it is
                % an object)
                value = R.content(idx).value;
                xplr.debuginfo('registry', 'key %s remove key (delete %s)', ...
                    fn_hash(R.content(idx).key,3), char(value))
                R.content(idx) = [];
                % -> and delete it! (if it is an object)
                if isobject(value), delete(value), end
                rm = true;
            else
                rm = false;
            end
            if nargout, removed = rm; end
            % Raise 'Empty' event if registry becomes empty
            if R.n == 0
                notify(R,'Empty')
            end
        end
        function value = getValue(R,key,newuser)
            value = []; % default output if key is not found
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = getValue(R,key{1});
                    if isempty(R1)
                        return
                    elseif ~isa(R1,'xplr.bankRegistry')
                        error 'key is invalid: encoutered a leave instead of a sub-registry'
                    end
                    value = getValue(R1,key(2:end),newuser);
                    return
                end
            end
            % Get value
            idx = getIndex(R,key);
            if isempty(idx), return, end
            value = R.content(idx).value;
            % No user: valid only with internal use, for nested
            % bankRegistry structure
            if nargin<3
                if ~isa(value,'xplr.bankRegistry'), error programming, end
                return
            end
            % Add new user
            idxuser = fn_find(newuser,R.content(idx).users(1,:),'first');
            if isempty(idxuser)
                xplr.debuginfo('registry', 'key %s add user %s', ...
                    fn_hash(R.content(idx).key,3), char(newuser))
                if isobject(newuser)
                    hld = addlistener(newuser,'ObjectBeingDestroyed', ...
                        @(u,e)R.unregister(key,newuser));
                else
                    hld = [];
                end
                R.content(idx).users(:,end+1) = {newuser; hld};
            else
                disp 'new user to register was already found in users list'
            end
        end
        function clear(R)
            R.content(:) = [];
        end
    end
    
end