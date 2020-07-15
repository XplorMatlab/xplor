classdef GraphNode < matlab.mixin.SetGet
    
    properties (Transient) %(Access='private')
        % objects the graph node is listening to
        listening = struct('object', {}, 'listener', {});
    end
    
    properties (Transient)
        id_graph_node
    end
    
    events
        test_event
    end
    
    % Constructor, display
    methods
        function self = GraphNode()
            self.id_graph_node = rand();
        end
        function delete(self)
            xplr.debug_info('GraphNode', ['delete ', class(self), num2str(floor(self.id_graph_node*1000), '%.3i')]);
            % when object is being deleted, make sure no more listener can
            % trigger actions on it
            for i = 1:length(self.listening)
                s = self.listening(i);
                brick.delete_valid(s.listener)
            end
        end
        function str = char(self)
            if isvalid(self)
                str = [class(self), num2str(floor(self.id_graph_node*1000), '%.3i')];
            else
                str = ['deleted ', class(self)];
            end
        end
    end
    
    % Listeners
    methods
        function add_listener(self, other, varargin)
            % function add_listener(self,other,add_listener arguments...)
            % function add_listener(self, other, listener)
            %---
            % Adding a listener on object other by using self.add_listener
            % ensures that when self will be deleted, the listener also
            % will be deleted and therefore not trigger actions on an
            % invalid object.
            if length(varargin) == 1
                % listener given in arguments
                listener = varargin{1};
            else
                % create listener here
                listener = addlistener(other, varargin{:});
            end
            self.listening(end + 1) = struct('object', other, 'listener', listener);
        end
        function add_listener_exclusive_pair(self, other, event1, callback1, event2, callback2)
            % function add_listener_exclusive_pair(self, other, event1, callback1, event2, callback2)
            %---
            % create two listeners inhibited when called by ach other to avoid loop
            % made to synchronize items with one central item
            % 
            % Inputs self_event and other_event can be event names, but also
            % 2-elements cell-array such as {'PostSet', PropertyName}.
            
            % listener self -> other (we cannot fully define it yet because
            % the second listener that will need to be inhibited is not
            % defined yet)
            if ~iscell(event1), event1 = {event1}; end
            listener1 = addlistener(self, event1{:}, @disp);
            other.listening(end+1) = struct('object', self, 'listener', listener1);
            % listener other -> self
            if ~iscell(event2), event2 = {event2}; end
            listener2 = addlistener(other, event2{:}, @(u,e)do_one_way_callback(callback2, u, e, listener1));
            self.listening(end+1) = struct('object', other, 'listener', listener2);
            % now we can define the callback of the first listener
            listener1.Callback = @(u,e)do_one_way_callback(callback1, u, e, listener2);
        end
        function disconnect(self, other)
            if ~isscalar(other)
                for obj = other
                    self.disconnect(obj)
                end
                return
            end
            % scan list of objects that self is listening, remove items
            % where self listens to other
            n = length(self.listening);
            rm = false(1, n);
            for i = 1:n
                s = self.listening(i);
                if s.object == other
                    rm(i) = true;
                    brick.delete_valid(s.listener)
                end
            end
            self.listening(rm) = [];
            % scan list of objects that other is listening, remove items
            % where other listens to self
            if ~isa(other, 'xplr.GraphNode') || ~isvalid(other), return, end
            n = length(other.listening);
            rm = false(1, n);
            for i = 1:n
                s = other.listening(i);
                if s.object == self
                    rm(i) = true;
                    brick.delete_valid(s.listener)
                end
            end
            other.listening(rm) = [];
        end
        function activate_connection(self, other, val)
            % Enable or disable a connection
            if ~isscalar(other)
                for obj = other
                    self.disconnect(obj)
                end
            end
            % scan list of objects that self is listening
            n = length(self.listening);
            for i = 1:n
                s = self.listening(i);
                if s.object == other
                    s.listener = val;
                end
            end
            % scan list of objects that other is listening
            if ~isa(other, 'xplr.GraphNode'), return, end
            n = length(other.listening);
            for i = 1:n
                s = other.listening(i);
                if s.object == self
                    s.listener = val;
                end
            end
        end
        function c = disable_connection(self, other)
            % Disable TEMPORARILY a connection: when c will be deleted, the
            % connection will automatically be re-enabled
            hl = [];
            % scan list of objects that self is listening
            n = length(self.listening);
            for i = 1:n
                s = self.listening(i);
                if s.object == other
                    hl(end+1) = s.listener;
                end
            end
            % scan list of objects that other is listening
            if isa(other, 'xplr.GraphNode')
                n = length(other.listening);
                for i = 1:n
                    s = other.listening(i);
                    if s.object == self
                        hl(end+1) = s.listener;
                    end
                end
            end
               
            % use disable_listener function (in brick)
            c = brick.disable_listener(hl);
        end
    end
    
    % Composition (i.e. components that exist iff object exists)
    methods
        function component = add_component(self, component)
            addlistener(self, 'ObjectBeingDestroyed', @(u,e)delete(component));
            addlistener(component, 'ObjectBeingDestroyed', @(u,e)delete(self));
            if nargout==0, clear component, end
        end
    end
     
    % Testing area
    methods
        function test(self)
            disp('hello')
            disp(self.testVar)
            notify(self, 'TestEvent')
        end
    end
    
end

function do_one_way_callback(callback, u, e, wayback_listener)
% Execute a callback, but first inhibit temporarily a given listener
    c = brick.disable_listener(wayback_listener);
    callback(u,e)
    % waybacklistener will be enable when c destroyed
end
