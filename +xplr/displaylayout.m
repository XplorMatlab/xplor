classdef displaylayout
    % function L = displaylayout(D)
    %---
    % This is a simple class that defines where each dimension is displayed
    % in XPLOR display.
    
    properties
        % parent display
        D
        % layout: the properties below are dimension identifiers (not
        % numbers)
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
            datadimID = [D.slice.header.dimID];
            if D.nd == 3 && doimage
                L.x = datadimID(1);
                L.y = datadimID(2);
                L.xy = datadimID(3);
            else
                L.x = datadimID(1:2:nd);
                L.y = datadimID(2:2:nd);
            end
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
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            F = {'x', 'y', 'ystatic', 'xy', 'yx'};
            for i=1:length(F)
                f = F{i};
                dimID = L.(f);
                [ok, dim] = ismember(dimID,datadimID);
                pos(dim(ok)) = {f}; 
            end
        end
        function pos = dim_location(L, dim)
            if ~isscalar(dim), error 'argument ''dim'' must be scalar', end
            dimID = L.D.slice.dimensionID(dim);
            F = {'x', 'y', 'ystatic', 'xy', 'yx'};
            for i=1:length(F)
                f = F{i};
                if any(L.(f) == dimID)
                    pos = f;
                    return
                end
            end
            error('no location found for dimension %g', dimID)
        end
        function s = dimensionNumber(L)
            % function s = dimensionNumber(L)
            %---
            % Return a structure with same fields x, y, ystatic, xy and yx
            % as L, but where values will be dimension numbers instead of
            % identifiers.
            s = struct;
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            F = {'x', 'y', 'ystatic', 'xy', 'yx'};
            for i=1:length(F)
                f = F{i};
                dimID = L.(f);
                [ok, dim] = ismember(dimID,datadimID);
                if ~all(ok)
                    error 'some dimension is not in data'
                end
                s.(f) = dim;
            end
        end
        function L2 = currentlayout(L)
            % function L2 = currentlayout(L)
            %---
            % keep from layout L only non-singleton dimensions that are
            % actually in L.D.slice
            L2 = L; % copy (displaylayout is not a handle class)
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            F = {'x' 'y' 'ystatic' 'xy' 'yx'};
            for i = 1:length(F)
                f = F{i};
                dimIDf = L.(f);      % dimension identifiers
                [present, dim] = ismember(dimIDf, datadimID);
                dim_present = dim(present); % dimension numbers of those actually present in the data
                singleton = ([datahead(dim_present).n]==1);
                L2.(f) = datadimID(dim_present(~singleton));
            end
        end
        function [L, L2] = updateLayout(L)
            % function [L, L2] = updateLayout(L)
            %---
            % Update layout upon slice change.
            % Keep locations of dimensions already present in L, 
            % remove dimensions that are not in L.D.slice anymore,
            % use some heuristic to choose locations of new dimensions
            
            L2 = L; % copy (displaylayout is not a handle class)
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            
            % dimensions already positionned + remove dimensions that
            % disappeared
            dimPositionned = false(1,L.D.nd);
            F = {'x' 'y' 'ystatic' 'xy' 'yx'};
            for i = 1:length(F)
                f = F{i};
                dimIDf = L.(f);      % dimension identifiers
                % mark among data dimensions those that are positionned
                dimPositionned = dimPositionned | ismember(datadimID, dimIDf);
                % and keep among L.(f) and L2.(f) only dimensions that are
                % still present in the data
                [present, dim] = ismember(dimIDf, datadimID);
                dim = dim(present); % dimension numbers of those actually present in the data
                dim_present = dim(present); % dimension numbers of those actually present in the data
                singleton = ([datahead(dim_present).n]==1);
                L.(f) = datadimID(dim_present);
                L2.(f) = datadimID(dim_present(~singleton));
            end
            
            % layout without the singleton dimensions
            L2 = L.currentlayout();
            
            % position missing dimensions
            for d = find(~dimPositionned)
                if datahead(d).n > 1
                    % insert in both L (which gathers all dimensions) and
                    % L2 (which keeps only displayed dimensions); try to
                    % keep L2 balanced between number of dimensions in x
                    % and y
                    if length(L2.y) <= length(L2.x)
                        L.y(end+1) = datadimID(d);
                        L2.y(end+1) = datadimID(d);
                    else
                        L.x(end+1) = datadimID(d);
                        L2.x(end+1) = datadimID(d);
                    end
                else
                    % insert only in L because dimension is not displayed;
                    % try to keep L balanced
                    if length(L.y) <= length(L.x)
                        L.y(end+1) = datadimID(d);
                    else
                        L.x(end+1) = datadimID(d);
                    end
                end
            end
        end
    end
    
end