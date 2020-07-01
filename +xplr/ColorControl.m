classdef ColorControl < xplr.GraphNode
   % color control
   
    properties
        hu      %
        values  %
        color   %
    end
    
    methods
        function u = ColorControl(varargin)
            u.values = {'white', 'red', 'blue'};
            u.hu = uicontrol('style', 'popupmenu', 'string', [{''}, u.values], ...
                'callback', @(hu, e)chg_color(u), ...
                varargin{:});
            u.color = [1, 1, 1];
        end
        function chg_color(u)
            val = get(u.hu, 'value');
            if val == 1, return, else val = val - 1; end
            set(u.hu, 'backgroundcolor', u.values{val}, 'value', 1)
            u.color = get(u.hu, 'backgroundcolor');
        end
    end
    
end
