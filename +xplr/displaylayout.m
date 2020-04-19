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
        mergeddata
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
            F = {'x', 'y', 'mergeddata', 'xy', 'yx'};
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
            F = {'x', 'y', 'mergeddata', 'xy', 'yx'};
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
            % Return a structure with same fields x, y, mergeddata, xy and yx
            % as L, but where values will be dimension numbers instead of
            % identifiers.
            s = struct;
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            F = {'x', 'y', 'mergeddata', 'xy', 'yx'};
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
            % keep from layout L only dimensions that are non-singleton
            L2 = L; % copy (displaylayout is not a handle class)
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            F = {'x' 'y' 'mergeddata' 'xy' 'yx'};
            for i = 1:length(F)
                f = F{i};
                dimIDf = L.(f);      % dimension identifiers
                [present, dimf] = ismember(dimIDf, datadimID);
                if ~all(present)
                    error 'programming: some layout dimensions are not present in the slice'
                end
                keep = ([datahead(dimf).n]>1);
                L2.(f) = datadimID(dimf(keep));
            end
        end
        function [L, L2] = updateLayout(L)
            % function [L, L2] = updateLayout(L)
            %---
            % Update layout upon slice change.
            % Keep locations of dimensions already present in L, 
            % remove dimensions that are not in L.D.slice anymore,
            % use some heuristic to choose locations of new dimensions
            % 
            % Input:
            % - L   layout
            %
            % Output:
            % - L   new layout
            % - L2  new layout, with singleton dimensions removed
            
            L2 = L; % copy (displaylayout is not a handle class)
            datahead = L.D.slice.header;
            datadimID = [datahead.dimID];
            
            % dimensions already positionned + remove dimensions that
            % disappeared
            dimPositionned = false(1,L.D.nd);
            F = {'x' 'y' 'mergeddata' 'xy' 'yx'};
            for i = 1:length(F)
                f = F{i};
                dimIDf = L.(f);      % dimension identifiers
                % mark among data dimensions those that are positionned
                dimPositionned = dimPositionned | ismember(datadimID, dimIDf);
                % and keep among L.(f) and L2.(f) only dimensions that are
                % still present in the data
                [present, dim] = ismember(dimIDf, datadimID);
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
        function [L, L2] = set_dim_location(L,dimID,location)
            % function [L, L2] = set_dim_location(L,dimID,location)
            %---
            % set new location of specific dimensions; locations of other
            % dimensions will automatically be adjusted
            %
            % Input:
            % - dimID       dimension identifier(s)
            % - location    character arrays: 'x', 'y', 'xy', 'yx' or
            %               'mergeddata', with, for 'x' and 'y' optional
            %               indication of index where to insert the
            %               location ('x0' = at the beginning = default,
            %               'x1' = next, 'x-1' = last)
            %               or cell array thereof for multiple dimensions
            %
            % Output:
            % - L   new layout
            % - L2  new layout, with singleton dimensions removed
            
            % input
            dimID = L.D.slice.dimensionID(dimID);
            if ~iscell(location), location = {location}; end
            if length(dimID) ~= length(location), error 'input lengths mismatch', end
            
            % remove these dimensions from the layout
            F = {'x' 'y' 'mergeddata' 'xy' 'yx'};
            for i = 1:length(F)
                f = F{i};
                L.(f)(ismember(L.(f), dimID)) = [];
            end
            
            % add them at the requested locations
            dimIDmove = []; 
            for i = 1:length(dimID)
                f = location{i};
                if ismember(f,{'xy' 'yx'})
                    % dimension that is already at 'xy' or 'yx' location if
                    % any will need to be moved somewhere else
                    dimIDmove = [L.xy L.yx]; 
                    [L.xy, L.yx] = deal([]);
                    % assign requested dimension to the requested location
                    L.(f) = dimID(i);
                else
                    % insert in either location 'x', 'y' or 'mergeddata' at a
                    % specific position
                    [f, index] = fn_regexptokens(f,'^(x|y|mergeddata)([\-0-9]*)$');
                    if isempty(index)
                        index = 0;
                    else
                        index = str2double(index);
                        if index<0
                            % negative index: count from the end
                            index = length(L.(f)) + 1 - index;
                        end
                    end
                    L.(f) = [L.(f)(1:index) dimID(i) L.(f)(index+1:end)];
                end
            end
            
            % layout without the singleton dimensions
            L2 = L.currentlayout();
            
            % dimension that needs to be moved somewhere else
            for d = dimIDmove
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