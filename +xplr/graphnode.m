classdef graphnode < matlab.mixin.SetGet
    
    properties %(Access='private')
        listening = struct('object',{},'listener',{});
    end
    
    properties
        idGraphNode
    end
    
    
    events
        TestEvent
    end
        
    
    methods
        function self = graphnode()
            self.idGraphNode = rand();
        end
        function delete(self)
            xplr.debuginfo('graphnode', 'delete %s', self)
            % when object is being deleted, make sure no more listener can
            % trigger actions on it
            for i = 1:length(self.listening)
                s = self.listening(i);
                deleteValid(s.listener)                
            end
        end
        function str = char(self)
            str = [class(self) num2str(floor(self.idGraphNode*1000),'%.3i')];
        end
        function addListener(self,other,varargin)
            % function addListener(self,other,addlistener arguments...)
            % function addListener(self,other,listener)
            if length(varargin) == 1
                % listener given in arguments
                listener = varargin{1};
            else
                % create listener here
                listener = addlistener(other,varargin{:});
            end
            self.listening(end+1) = struct('object',other,'listener',listener);
        end
        function addListenerExclusivePair(self,other,arg1,callback1,arg2,callback2)
            % create two listeners inhibited when called by each other to avoid loop
            % made to synchronize items with one central item
            
            % listener self -> other (we cannot fully define it yet because
            % the second listener that will need to be inhibited is not
            % defined yet)
            if ~iscell(arg1), arg1 = {arg1}; end
            listener1 = addlistener(self,arg1{:},@disp);
            other.listening(end+1) = struct('object',self,'listener',listener1);
            % listener other -> self
            if ~iscell(arg2), arg2 = {arg2}; end
            listener2 = addlistener(other,arg2{:},@(u,e)doOneWayCallback(callback2,u,e,listener1));
            self.listening(end+1) = struct('object',other,'listener',listener2);
            % now we can define the callback of the first listener
            listener1.Callback = @(u,e)doOneWayCallback(callback1,u,e,listener2);
        end
        function disconnect(self,other)
            if ~isscalar(other)
                for obj = other
                    self.disconnect(obj)
                end
                return
            end
            % scan list of objects that self is listening, remove items
            % where self listens to other
            n = length(self.listening);
            rm = false(1,n);
            for i = 1:n
                s = self.listening(i);
                if s.object == other
                    rm(i) = true;
                    deleteValid(s.listener)
                end
            end
            self.listening(rm) = [];
            % scan list of objects that other is listening, remove items
            % where other listens to self
            if ~isa(other,'xplr.graphnode'), return, end
            n = length(other.listening);
            rm = false(1,n);
            for i = 1:n
                s = other.listening(i);
                if s.object == self
                    rm(i) = true;
                    deleteValid(s.listener)
                end
            end
            other.listening(rm) = [];
        end
        function activateConnection(self,other,val)
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
                    enableListener(s.listener, val)
                end
            end
            % scan list of objects that other is listening
            if ~isa(other,'xplr.graphnode'), return, end
            n = length(other.listening);
            for i = 1:n
                s = other.listening(i);
                if s.object == self
                    enableListener(s.listener, val)
                end
            end
        end
        function c = disableConnection(self,other)
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
            if isa(other,'xplr.graphnode')
                n = length(other.listening);
                for i = 1:n
                    s = other.listening(i);
                    if s.object == self
                        hl(end+1) = s.listener;
                    end
                end
            end
               
            % use disableListener function (in brick)
            c = disableListener(hl);
        end
    end
     
    % Testing area
    methods
        function test(self)
            disp('hello')
            disp(self.testVar)
            notify(self,'TestEvent')
        end
    end
    
end

function doOneWayCallback(callback,u,e,waybacklistener)
% Execute a callback, but first inhibit temporarily a given listener
    c = disableListener(waybacklistener);
    callback(u,e)
    % waybacklistener will be enable when c destroyed
end