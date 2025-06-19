classdef PointListener < xplr.GraphNode
% Creates and returns a listner that executes callback when current point
% coordinates in specified dimension(s) change
%---
% function listener = point_listener(obj, dim, callback[, world])
%---
% Input:
% - obj         a View, its Display or Navigation object
% - dim         dimension(s) to watch for
% - callback    function to execute with prototype fun(coord)
% - world       whether input coordonates are world coordonates or data
%               coordinates [default: false, i.e. use data coordinates]


properties
    P
    listener
    callback
end


methods
    function L = PointListener(obj, dim, callback, world)      
        if nargin>3 && world
            error 'not implemented yet'
        end
        
        % Get Navigation object
        if isa(obj, 'xplr.DisplayNavigation')
            N = obj;
        elseif isa(obj, 'xplr.ViewDisplay')
            N = obj.navigation;
        elseif isa(obj, 'xplr.View')
            N = obj.D.navigation;
        else
            error 'object must be a View, Display or Navigation object'
        end
        
        % Find point filter for selected dimension
        d = N.D.slice.dimension_number(dim); % convert dimension to number
        L.P = N.point_filters{d};
        % register L as a new user of this point filter, to make sure the
        % point filter won't be deleted when some other xplor objects will
        % be deleted
        L.P.add_user(L)
        
        % Memorize callback and try to execute it once
        L.callback = callback;
        msg = 'Try to execute callback once.';
        disp(msg)
        try
            L.callback(L.P.index)
        catch
            fprintf(repmat('\b', 1, length(msg)+1));
            disp('Callback failed. Not creating listener.')
            % repeat the error so that it can be debugged
            L.callback(L.P.index)
        end
            
        % Add listener
        L.add_listener(L.P, 'changed_operation', @(u,e)moved_point(L, e));
        fprintf(repmat('\b', 1, length(msg)+1));
        disp(['Listener created to changes in ''' L.P.header_in.label ...
            ''' dimension. Callback executed once.'])
        
        % No output
        if nargout == 0
            clear L
        end
    end
    
    function moved_point(L, e)     
        if ~e.chg_ij
            % slice coordinate did not change
            return
        end
        
        % Try execute callback
        try
            L.callback(L.P.index)
        catch err
            % If failure, maybe some object necessary to callback execution
            % was deleted. Analyze error message to determine whether it is
            % the case, in such case silently delete the listner.
            msg_start = ['Callback to change in ''' L.P.header_in.label ...
                ''' dimension failed'];
            if strcmp(err.message, 'Invalid or deleted object.')
                delete(L)
                disp([msg_start ' because an object became invalid. Removed listener.'])
            else
                answer = questdlg( ...
                    {'An error occured while executing point listener callback:' ...
                    err.message ...
                    'Do you want to delete the listnener?'});
                if strcmp(answer, 'Yes')
                    delete(L)
                else
                    % Execute callback again for user to see the error, maybe
                    % debug it.
                    L.callback(L.P.index)
                end
            end
        end     
    end
end

end