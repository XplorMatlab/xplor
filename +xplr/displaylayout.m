classdef displaylayout
    % function L = displaylayout(D)
    %---
    % This is a simple class that defines where each dimension is displayed
    % in XPLOR display.
    
    properties
        % parent display
        D
        % layout 
        x
        y
        ystatic
        xy
        yx
    end
    properties (Dependent, SetAccess='private')
        dim_locations
        griddisplay
    end
    
    % Create
    methods
        function L = displaylayout(D)
            % parent display
            L.D = D;
            
            % choose a default organization for display D
            nd = D.nd;
            doimage = strcmp(D.displaymode,'image');
            if D.nd == 3 && doimage
                L.x = 1;
                L.y = 2;
                L.xy = 3;
            else
                L.x = 1:2:nd;
                L.y = 2:2:nd;
            end
        end
    end
    methods (Static)
        function L = disconnectedLayout(sz,displaymode)
            % emulate a viewdisplay object D such that all displaylayout
            % methods below will work
            D = struct;
            D.nd = length(sz);
            D.slice.sz = sz;
            D.displaymode = displaymode;
            L = xplr.displaylayout(D);
        end
    end
    
    % Simple get dependent
    methods
        function d = get.griddisplay(L)
            d = [L.xy L.yx];
        end
    end
    
    % Manipulation
    methods
        function pos = get.dim_locations(L)
            pos = cell(1,L.D.nd);
            F = {'x', 'y', 'ystatic', 'xy', 'yx'};
            for i=1:length(F)
                f = F{i};
                pos(L.(f)) = {f}; 
            end
        end
        function L2 = strip_singleton(L, D)
            L2 = L; % copy (displaylayout is not a handle class)
            singleton = (L.D.slice.sz == 1);
            L2.x(singleton(L2.x)) = [];
            L2.y(singleton(L2.y)) = [];
            L2.ystatic(singleton(L2.ystatic)) = [];
            L2.xy(singleton(L2.xy)) = [];
            L2.yx(singleton(L2.yx)) = [];
        end
    end
    
end