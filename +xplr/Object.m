classdef Object < matlab.mixin.SetGet
    
    properties (Transient)
        object_id
    end
    
    % Constructor, display
    methods
        function self = Object()
            self.object_id = rand();
            xplr.debug_info('Object', ['create ', char(self)]);
        end
        function delete(self)
            xplr.debug_info('Object', ['delete ', char(self)]);
        end
        function str = char(self)
            cls = class(self);
            cls = cls(6:end); % remove 'xplr.'
            try
                str = [cls, num2str(floor(self.object_id*1000), '%.3i')];
            catch
                str = ['deleted ', cls];
            end
        end
    end
    
end

