classdef BankRegistry < handle
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
        content = struct('key', cell(1, 0), 'value', cell(1,0), 'users', cell(1, 0));
    end
    
    events
        Empty
    end
    
    methods
        function n = get.n(R)
            n = length(R.content);
        end
        function str = char(R)
            str = 'BankRegistry';
        end
    end
    
    methods (Access='private')
        function idx = get_index(R, key)
            idx = fn_find(key, {R.content.key}, 'first');
        end
    end
    methods
        function register(R, key, value, user)
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = get_value(R, key{1});
                    if isempty(R1)
                        % create sub-register
                        R1 = xplr.BankRegistry();
                        idx = R.n + 1;
                        R.content(idx) = struct('key', key{1}, 'value', R1, 'users', {cell(1, 0)});
                        % watch for it to become empty
                        addlistener(R1, 'Empty', @(u,e)unregister(R,key{1}));
                    elseif ~isa(R1, 'xplr.BankRegistry')
                        error 'key is invalid: encountered a leave instead of a sub-registry'
                    end
                    register(R1, key(2:end), value, user)
                    return
                end
            end
            % Verify that no entry exists yet
            if ~isempty(get_index(R, key))
                if ischar(key) || isnumeric(key)
                    error('Register already has a value for key ''%s''.', num2str(key))
                else
                    error('Register already has a value for the specified key.')
                end
            end
            % Create new entry
            idx = R.n+1;
            if isempty(user)
                users = cell(2, 0);
            else
                if isobject(user)
                    hld = addlistener(user, 'ObjectBeingDestroyed', @(u,e)R.unregister(key, user));
                else
                    hld = [];
                end
                users = {user; hld};
            end
            R.content(idx) = struct('key', key, 'value', value, 'users', {users});
            xplr.debug_info('registry', 'key %s 1st user %s (value %s)', ...
                fn_hash(key, 3), char(user), char(value))
        end
        function removed = unregister(R, key, user)
            if ~isvalid(R), return, end % R has been deleted already
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = get_value(R, key{1});
                    if isempty(R1)
                        disp 'sub-registry where to unregister was not found'
                        if nargout, removed = false; end
                        return
                    elseif ~isa(R1, 'xplr.BankRegistry')
                        error 'key is invalid: encoutered a leave instead of a sub-registry'
                    end
                    if nargin<3
                        rm = unregister(R1, key(2:end)); % note that if R1 becomes empty, it will emit 'Empty' and automatically get removed from R and be deleted
                    else
                        rm = unregister(R1, key(2:end), user); % note that if R1 becomes empty, it will emit 'Empty' and automatically get removed from R and be deleted
                    end
                    if nargout, removed = rm; end
                    return
                end
            end
            % Value index
            idx = get_index(R, key);
            if isempty(idx)
                xplr.debug_info('registry', 'key %s was not found', fn_hash(key,3))
                disp 'entry to unregister was not found'
                if nargout, removed = false; end
                return
            end
            % Unregister value if no more users
            if nargin<3
                % remove value regardless of there are users or not
                R.content(idx).users = [];
            else
                idxuser = fn_find(user, R.content(idx).users(1, :), 'first');
                if isempty(idxuser)
                    disp 'user to unregister was not found in users list'
                else
                    c = R.content(idx);
                    xplr.debug_info('registry', 'key %s rmv user %s', ...
                        fn_hash(c.key, 3), char(c.users{1, idxuser}))
                    hld = R.content(idx).users{2, idxuser}; % handle of listener listening for object deletion
                    delete(hld)
                    R.content(idx).users(:, idxuser) = [];
                end
            end
            if isempty(R.content(idx).users)
                % no more user -> unregister value and delete if (if it is
                % an object)
                value = R.content(idx).value;
                xplr.debug_info('registry', 'key %s remove key (delete %s)', ...
                    fn_hash(R.content(idx).key, 3), char(value))
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
                notify(R, 'Empty')
            end
        end
        function value = get_value(R, key, new_user)
            value = []; % default output if key is not found
            % Multi-level?
            if iscell(key)
                if isscalar(key)
                    key = key{1};
                else
                    R1 = get_value(R, key{1});
                    if isempty(R1)
                        return
                    elseif ~isa(R1, 'xplr.BankRegistry')
                        error 'key is invalid: encoutered a leave instead of a sub-registry'
                    end
                    value = get_value(R1, key(2:end), new_user);
                    return
                end
            end
            % Get value
            idx = get_index(R, key);
            if isempty(idx), return, end
            value = R.content(idx).value;
            % No user, not even empty one: valid only with internal use,
            % for nested bankRegistry structure
            if nargin<3
                if ~isa(value, 'xplr.BankRegistry'), error programming, end
                return
            end
            % No user to register: this is a dangerous situation as the
            % registry entry and the value can be deleted at any time, upon
            % unregistration of another user
            if isempty(new_user)
                xplr.debug_info('registry', 'key %s accessed without user',fn_hash(R.content(idx).key,3))
                return
            end
            % Add new user
            idx_user = fn_find(new_user, R.content(idx).users(1, :), 'first');
            if isempty(idx_user)
                xplr.debug_info('registry', 'key %s add user %s', ...
                    fn_hash(R.content(idx).key, 3), char(new_user))
                if isobject(new_user)
                    hld = addlistener(new_user,'ObjectBeingDestroyed', ...
                        @(u,e)R.unregister(key,new_user));
                else
                    hld = [];
                end
                R.content(idx).users(:, end+1) = {new_user; hld};
            else
                disp 'new user to register was already found in users list'
            end
        end
        function clear(R)
            R.content(:) = [];
        end
    end
    
end
