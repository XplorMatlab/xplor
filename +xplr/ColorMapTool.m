classdef ColorMapTool < xplr.GraphNode
    
    properties (SetObservable=true, AbortSet)
        c_map_def
        do_nonlinear = false
        invert_map = false
    end
    properties (SetObservable=true, AbortSet, SetAccess='private')
        c_map
    end
    properties (Access='private')
        menu
        nonlinear_fun_editor
    end
    
    events
        ChangedColorMap
    end
   
    % Constructor, destructor, menu
    methods
        function C = ColorMapTool(D)
            % Build menu
            buildmenu(C, D);
            % Set to default: jet
            C.c_map_def = 'jet';
        end
        function delete(C)
            delete@xplr.GraphNode(C)
            if ~isprop(C, 'menu'), return, end
            brick.delete_valid(C.menu)
        end
        function buildmenu(C, D)
            hf = D.V.hf;
            
            % Create menu
            if isempty(C.menu)
                m = uimenu('parent', hf, 'label', 'Color Map');
                C.menu = m;
            else
                m = C.menu;
                if nargin >= 2 && get(m, 'parent') ~= hf, error 'existing menu does not have the same parent figure', end
                delete(get(m, 'children'))
            end
            
            % List of possible color maps
            map_names = brick.row(fieldnames(xplr.ColorMapTool.available_color_maps()));
            brick.propcontrol(C, 'c_map_def', ['menugroup', map_names, 'user...'], m);
            
            % Apply non-linear function to values before coloring
            brick.propcontrol(C, 'invert_map', 'menu', ...
                {'parent', m, 'label', 'Invert map', 'separator', 'on'});
            brick.propcontrol(C, 'do_nonlinear', 'menu', ...
                {'parent', m, 'label', 'Apply nonlinear function before coloring'});
            
            % Control visibility depending on dislay mode
            set(C.menu, 'visible', brick.switch_case(D.display_mode, 'image', 'on', 'off'));
            brick.connect_listener(D, C, 'display_mode_', 'PostSet', ...
                @(u,e)set(C.menu, 'visible', brick.switch_case(D.display_mode, 'image', 'on', 'off')));
        end
    end
    
    % Available color maps
    methods (Static)
        function maps = available_color_maps()
            maps = struct( ...
                'gray',     @gray, ...
                'jet',      @jet, ...
                'parula',   @parula, ...
                'hot',      @hot, ...
                'mapgeog',  @colormaps.mapgeog, ...
                'mapgeogclip',  @colormaps.mapgeogclip, ...
                'mapclip',  @colormaps.mapclip, ...
                'mapcliphigh',  @colormaps.mapcliphigh, ...
                'mapcliplow',   @colormaps.mapcliplow, ...
                'red',      @colormaps.red, ...
                'green',    @colormaps.green, ...
                'blue',     @colormaps.blue, ...
                'black_red',    @colormaps.black_red, ...
                'black_green',  @colormaps.black_green, ...
                'black_blue',   @colormaps.black_blue, ...
                'white_red',    @colormaps.white_red, ...
                'white_green',  @colormaps.white_green, ...
                'white_blue',   @colormaps.white_blue, ...
                'white_black',  @colormaps.white_black, ...
                'maporient',    @colormaps.maporient ...
                );
        end
    end
    
    % Color map and color map name: be careful with infinite loops!
    methods
        function set.c_map_def(C, x)
            % Handle colormap names and 'user..' command
            if ischar(x)
                name = x;
                switch name
                    case 'user...'
                        answer = inputdlg('Define color map here:', '', 1, {'rand(256,3)'});
                        if isempty(answer), return, end
                        try
                            x = evalin('base', answer{1});
                        catch
                            errordlg('Could not evaluate command to a valid color map', '', 'modal')
                            return
                        end
                    otherwise
                        maps = xplr.ColorMapTool.available_color_maps();
                        if isfield(maps, name)
                            x = feval(maps.(name), 256);
                        else
                            errordlg(sprintf('''%s'' is not the name of a valid color map', name), '', 'modal')
                            return
                        end
                end
            else
                name = 'user';
            end
            % Set color map 
            if ~(isnumeric(x) && ismatrix(x) && size(x, 2) == 3 && all(x(:) >= 0 & x(:) <= 1))
                error('not a valid color map definition')
            end
            C.c_map = x; %#ok<MCSUP>
            % Set name
            C.c_map_def = name;
            % Notify
            notify(C, 'ChangedColorMap')
        end
    end
    
    % Nonlinear function edition: we apply a nonlinear function before
    % applying the colormap, this allows for example in a black & white
    % image to enhance low contrast in the dark range, etc.
    methods
        function set.do_nonlinear(C, value)
            C.do_nonlinear = logical(value);
            % is we want to apply nonlinear function, we create an object
            % of class signaleditor to control the parameters of this
            % function
            if value && (isempty(C.nonlinear_fun_editor) || ~isvalid(C.nonlinear_fun_editor))
                C.nonlinear_fun_editor = brick.signaleditor([1, 256], [0, 1], ...
                    @(x)notify(C, 'ChangedColorMap'), 'monotonous', 'min', 0, 'max', 1);
                % if the editor is closed, we stop applying the nonlinear
                % function
                C.add_listener(C.nonlinear_fun_editor, 'ObjectBeingDestroyed', @(u,e)set(C, 'do_nonlinear', false))
            else
                notify(C, 'ChangedColorMap')
            end
        end
        function set.invert_map(C, value)
            C.invert_map = value;
            % update diplay
            notify(C, 'ChangedColorMap')
        end
    end
    
    % Apply colormap to image
    methods
        function im = color_image(C, xi, clipi)
            % vectorize image
            [nx, ny, nc] = size(xi);
            xi = reshape(xi, [nx*ny, nc]);
            idxnan = any(isnan(xi), 2);
            % clip data; note that clipi can have 3 clipping ranges for the
            % 3 image color channels
            clipi = brick.matrix(clipi);
            xi = brick.div(brick.subtract(xi, clipi(1,:)), diff(clipi));
            xi = max(0, min(1, xi));
            % apply nonlinear function if any
            if C.do_nonlinear
                xi = C.nonlinear_fun_editor.interp(1 + 255*xi);
            end
            % color the image and replace NaNs by a specific color
            if size(xi, 2) == 1
                xi(idxnan) = 0;
                if C.invert_map
                    xi = 256 - round(xi*(size(C.c_map, 1) - 1));
                else
                    xi = 1 + round(xi*(size(C.c_map, 1) - 1));
                end
                im = C.c_map(xi(:), :);
            else
                im = xi;
            end
            nanvalue = .95;
            im(idxnan, :) = nanvalue;
            % final reshape restores original image size
            im = reshape(im, [nx, ny, size(im, 2)]);
        end
    end
end
