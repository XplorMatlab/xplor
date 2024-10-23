classdef white2alpha < brick.interface
    %WHITE2ALPHA Convert image white background to transparency
    %---
    % function white2alpha(image)
    % function white2alpha([filename|'clipboard'])
    %---
    % Convert image white background to transparency
    % 
    % See also brick.show_image, brick.combine_images, brick.montage,
    % brick.image_anonymous

    properties (SetAccess='private')
        controls
        filename
        original_input
        bg_color
        input
        lumspecial
        saturation
        checker
        truecolor
        alpha
        image_subpart
        nx
        ny
        nc
    end
    properties (SetObservable)
        dosubimage = '';
    end
    properties (Access='private')
        menu_def_subpart
        im = struct;
        borders
        image_subpart_poly
        image_subpart_mark
    end

    methods
        function X = white2alpha(varargin)

            X = X@brick.interface(917,'WHITE2ALPHA');

            % Graphic objects
            X.grob = struct;
            X.grob.controls = uipanel;
            X.grob.input = axes;
            X.grob.masks = axes;
            X.grob.truecolor = axes;
            X.grob.alpha = axes;
            X.grob.result = axes;
            X.interface_end()

            % Controls
            s = struct(...
                'special__max__alpha__method', {false 'logical'}, ...
                'outside__max__luminance',  {.99    'slider .3 1 .01 %.2f < ~special__max__alpha__method'}, ...
                'holes',                    {false  'logical < ~special__max__alpha__method'}, ...
                'flat__colors',             {true   'logical < ~special__max__alpha__method'}, ...
                'border__typical__width',   {0.5    'logslider 0 3 .01 < ~special__max__alpha__method'}, ...
                'flat__color__tolerance',   {.01    'slider 0 .1 .005 %.2f < flat__colors ~special__max__alpha__method'}, ...
                'border__max__luminance',   {.5     'slider 0 1 .005 %.2f < ~flat__colors ~special__max__alpha__method'}, ...
                'border__max__saturation',  {.5     'slider 0 1 .005 %.2f < ~flat__colors ~special__max__alpha__method'}, ...
                'final__alpha__smooth',            {false  'logical < ~special__max__alpha__method ~special__max__alpha__method'}, ...
                'true__color__smoothing',   {0      'slider 0 1 < ~final__alpha__smooth ~special__max__alpha__method'});
            X.controls = brick.control(s,@(s)X.action(s),X.grob.controls);

            % Load image and perform conversion
            X.load_image(varargin{:})
        end
        function load_image(X, a, varargin)
            % Options
            do_sub_region = brick.flags({'subregion'}, varargin);
            
            % Image
            if nargin>=2 && strcmp(a, 'clipboard')
                % try pasting from clipboard; if there is no image in
                % clipboard, result will be empty
                a = permute(imclipboard('paste'), [2 1 3]);
                switch class(a)
                    case 'double'
                    case 'uint8'
                        % convert to double
                        a = double(a) / 255;
                    otherwise
                        error('type %s not handled yet, please edit code', class(a))
                end
            end
            if nargin<2 || isempty(a)
                a = brick.getfile('*','Select image');
                if isequal(a, 0), return, end
            end
            while ischar(a)
                X.filename = a;
                [a, alph] = brick.readimg(a);
                switch class(a)
                    case 'double'
                        ok = all(alph(:)==1);
                    case 'uint8'
                        % convert to double
                        a = double(a) / 255;
                        ok = all(alph(:)==uint8(255));
                        alph = double(alph) / 255;
                    otherwise
                        error('type %s not handled yet, please edit code', class(a))
                end
                if ~ok
                    answer = questdlg('Image already has an alpha channel, what do you want to do?', ...
                        'brick.white2alpha','Use image','Use image weighted by alpha','Other image','Use image weighted by alpha');
                    if isempty(answer)
                        % interrupt, close window
                        close(X.hf)
                        return
                    end
                    switch answer
                        case 'Use image'
                            % use a as is
                        case 'Use image weighted by alpha'
                            a = brick.mult(a, alph)  + (1-alph);
                        case 'Other image'
                            a = brick.getfile('*','Select image');
                    end
                end
            end
            X.original_input = a;
            X.bg_color = [1 1 1];
            X.input = a;
            
            % Select sub-region
            if do_sub_region
                select_sub_region(X)
            end
            
            % Some precomputations
            precomputations(X)
        end
        function precomputations(X)
            [X.nx, X.ny, X.nc] = size(X.input);
            if X.nc == 3
                [~, X.saturation, X.lumspecial] = rgb2hsv(X.input);
                % replace luminance by darkest channel!
                X.lumspecial = min(X.input,[],3);
            elseif X.nc == 1
                X.saturation = 0;
                X.lumspecial = X.input;
            else
                error 'input image must have 1 or 3 color channels'
            end
            step = round(mean([X.nx X.ny])/10);
            xstripes = mod(floor((0:X.nx-1)'/step),2);
            ystripes = mod(floor((0:X.ny-1)/step),2);
            X.checker = bsxfun(@xor,xstripes,ystripes);
            X.checker = 1 - (1 - X.checker)*1;

            % Show images
            X.show_images()

            % Perform conversion
            X.performconversion()
        end
        function show_images(X)
            colormap(X.hf,gray(256))
            X.im.raw = imagesc(permute(X.input,[2 1 3]),'parent',X.grob.input,[0 1]);
            axis(X.grob.input,'image')
            set(X.grob.input,'xtick',[],'ytick',[],'box','on')
            title(X.grob.input,'Input')
            X.im.masks = imagesc(permute(X.input,[2 1 3]),'parent',X.grob.masks);
            axis(X.grob.masks,'image')
            set(X.grob.masks,'xtick',[],'ytick',[],'box','on')
            title(X.grob.masks,'Outside & border masks')
            X.im.truecolor = imagesc(permute(X.input,[2 1 3]),'parent',X.grob.truecolor,[0 1]);
            axis(X.grob.truecolor,'image')
            set(X.grob.truecolor,'xtick',[],'ytick',[],'box','on')
            title(X.grob.truecolor,'True color')
            X.im.alpha = imagesc(permute(X.input,[2 1 3]),'parent',X.grob.alpha,[0 1]);
            axis(X.grob.alpha,'image')
            set(X.grob.alpha,'xtick',[],'ytick',[],'box','on')
            title(X.grob.alpha,'Alpha')
            X.im.result = imagesc(permute(X.input,[2 1 3]),'parent',X.grob.result,[0 1]);
            axis(X.grob.result,'image')
            set(X.grob.result,'xtick',[],'ytick',[],'box','on')
            title(X.grob.result,'Result')
            brick.imvalue('image')
        end
        function select_sub_region(X)
            a = X.input;
            % select mask
            mask = brick.maskselect(a, 'free');
            % all pixels outside mask are white
            [n_x, n_y, ~] = size(a);
            a = brick.imvect(a, 'vector');
            a(~mask, :) = 1;
            a = brick.imvect(a, [n_x n_y]);
            % crop image to mask sides
            a = a(any(mask, 2), any(mask, 1), :);
            % store result
            X.input = a;
        end
        function select_background_color(X)
            
            % select background color
            a = X.original_input;
            hf = brick.figure('Select background color', 'tag', 'no');
            imshow(permute(a,[2 1 3]));
            ok = false;
            while ~ok
                waitforbuttonpress()
                p = get(gca, 'currentpoint');
                p = round(p(1, 1:2));
                ok = all(p >= 1 & p <= [size(a,1) size(a,2)]);
            end
            X.bg_color = brick.row(a(p(1),p(2),:));
            close(hf)
            
            % adapt image for a white background 
            % -> make a rough estimation of alpha and true color and
            % reimprint these true colors on a white background
            [xalpha, xtruecolor] = X.max_alpha_estimation();
            X.input = brick.add(brick.mult(xalpha, xtruecolor), 1-xalpha);
            
            % perform precomputations and display
            X.precomputations()
        end
        function init_menus(X)
            init_menus@brick.interface(X)

            % Load image
            m = X.menus.interface;
            uimenu(m,'label','Load image...','separator','on', ...
                'callback',@(u,e)X.load_image())
            uimenu(m,'label','Load image (sub-region)...'   , ...
                'callback',@(u,e)X.load_image([], 'subregion'))
            uimenu(m,'label','Load image from clipboard', ...
                'callback',@(u,e)X.load_image('clipboard'))
            uimenu(m,'label','Select original background color', ...
                'callback',@(u,e)X.select_background_color())
            uimenu(m,'label','Save result to file...','separator','on', ...
                'callback',@(u,e)X.save())

            % Image sub-part
            m = uimenu(X.hf,'label','Sub-Image');
            X.menus.image_sub_part = m;
            brick.propcontrol(X,'dosubimage', ...
                {'menugroup' {'inside' 'outside' ''} {'inside selection' 'outside selection'}}, ...
                'parent',m,'label','Apply to image sub-part');
            X.menu_def_subpart = uimenu(m, ...
                'label','Define new selection', ...
                'enable', brick.onoff(~isempty(X.dosubimage)), ...
                'callback',@(u,e)set_image_subpart(X));
            uimenu(m,'label','Reset display','separator','on', ...
                'callback',@(u,e)X.reset_display())

        end
        function set.dosubimage(X,val)
            X.dosubimage = val;
            b = ~isempty(val);
            set(X.menu_def_subpart,'enable',brick.onoff(b))
            if b && isempty(X.image_subpart)
                set_image_subpart(X)
            end
            set(X.image_subpart_mark,'visible',brick.onoff(b))
        end
        function set_image_subpart(X)
            brick.delete_valid(X.image_subpart_mark)
            poly = brick.mouse(X.grob.result,'poly','select image sub-part');
            poly = poly(:,[1:end 1]);
            X.image_subpart_poly = poly;
            X.show_image_subpart()
            X.image_subpart = brick.poly2mask(poly(1,:),poly(2,:),X.nx,X.ny);
            % restore view limits if they were modified
            axis(X.grob.result, axis(X.grob.input))
        end
        function show_image_subpart(X)
            brick.delete_valid(X.image_subpart_mark)
            X.image_subpart_mark = brick.drawpoly(X.image_subpart_poly, ...
                'parent',X.grob.result,'color','w', ...
                'visible',brick.onoff(~isempty(X.dosubimage)));
        end
        function reset_display(X)
            set(X.hf,'WindowButtonMotionFcn','')
            X.show_images()
            X.display_intermediary_result()
            X.display_final_result()
            X.show_image_subpart()
        end
        function action(X, ~)
            % parameter change
            X.performconversion()
        end
        function [xalpha, xtruecolor] = max_alpha_estimation(X)
            % -> we note that a = alpha x + (1-alpha) bg
            % hence x = (a - (1-alpha) bg) / alpha
            % we must have x >= 0 and x <= 1
            % we obtain alpha >= (bg-a) / bg and alpha >= (a-bg) / (1-bg)
            a = X.original_input;
            a = brick.imvect(a);  % (nx*ny) x 3
            bg = repmat(X.bg_color, [X.nx*X.ny, 1]);
            alpha_min = zeros(X.nx*X.ny, 3);
            a2 = brick.subtract(a, bg);
            sub = (a2 < 0);
            alpha_min(sub) = -a2(sub) ./ bg(sub);
            sub = (a2 > 0);
            alpha_min(sub) = a2(sub) ./ (1 - bg(sub));
            xalpha = max(alpha_min, [], 2); % (nx*ny) vector
            xalpha3 = repmat(xalpha, [1 3]);
            xtruecolor = (a - (1-xalpha3).*bg) ./ xalpha3;
            xtruecolor(xalpha3==0) = 0;
            xalpha = brick.imvect(xalpha, [X.nx X.ny]);
            xtruecolor = brick.imvect(xtruecolor, [X.nx X.ny]);
        end
        function performconversion(X)
            c = brick.watch(X.hf);
            
            % Special: max alpha method
            % maximal alpha possible depends on the darker channel for each
            % pixel
            if X.controls.special__max__alpha__method
                [xalpha, xtruecolor] = X.max_alpha_estimation();
                % save xalpha and xtruecolor, possible only for sub-part of
                % image
                X.save_result(xalpha, xtruecolor);
                % in borders, mark pixels with zero alpha in yellow, and
                % other pixels with non-ones alpha in black
                a = brick.imvect(X.input,'vector');
                a(xalpha<1, :) = 0;
                a(xalpha==0, 1:2) = 1;
                a = brick.imvect(a,[X.nx X.ny],'image');
                X.borders = a;
                % display result
                X.display_intermediary_result('original')
                X.display_final_result()
                return
            end

            % Inside and outside masks
            outside = (X.lumspecial > X.controls.outside__max__luminance);
            if ~X.controls.holes
                % consider only connected components of "outside" that
                % contain some border of the image
                labels = bwlabel(outside);
                ok_labels = unique([brick.row(labels([1 end],:)) brick.row(labels(:,[1 end]))]);
                outside = outside & ismember(labels,ok_labels);
            end
            if X.controls.flat__colors
                darkness = 1 - X.lumspecial;
                spread = round(X.controls.border__typical__width);
                local_max_darkness = imdilate(darkness,true(spread));
                inside = (darkness > local_max_darkness * (1-X.controls.flat__color__tolerance));
                inside = inside & ~outside;
            else
                inside = (X.lumspecial < X.controls.border__max__luminance);
                if X.nc == 3
                    inside = inside | (X.saturation > X.controls.border__max__saturation);
                end
            end
            border = ~(inside | outside);

            % border must touch outside
            labels = bwlabel(border);
            ok_labels = unique(labels(bwmorph(outside,'dilate')));
            border = border & ismember(labels,ok_labels);

            % border must not be further than a given distance to outside
            border = border & bwmorph(outside,'dilate',round(X.controls.border__typical__width));
            inside = ~(border | outside);

            a = brick.imvect(X.input,'vector');
            a(outside(:),1:2) = 1;
            a(outside(:),3) = 0;
            a(border,:) = 0;
            a = brick.imvect(a,[X.nx X.ny],'image');
            X.borders = a;
            X.display_intermediary_result()

            % Stop here if slider is being moved
            if X.controls.sliderscrolling
                return
            end

            % True color of semitransparent pixels obtained by smoothing
            % of inside
            sigma = 4 * X.controls.border__typical__width;
            negative_inside = brick.mult(1 - X.input, inside);
            negative_inside_smooth = brick.filt(negative_inside,sigma,'l',[1 2]);
            inside_mask_smooth = brick.filt(inside,sigma,'l',[1 2]);
            negative_inside_spread = brick.div(negative_inside_smooth, inside_mask_smooth);
            inside_spread = brick.clip(1 - negative_inside_spread, [0 1]);
            inside_spread = brick.imvect(inside_spread,'vector');
            xtruecolor = brick.imvect(X.input,'vector');
            truecolor_border = inside_spread(border, :);

            % Transparency of border pixels obtained as the ratio of input
            % image darkness as compared to true color darkness
            input_border = brick.imvect(X.input,border,'vector');
            darkness = 1 - mean(input_border,2);
            truecolor_border_darkness = 1 - mean(truecolor_border,2);
            alpha_border = brick.clip(darkness ./ truecolor_border_darkness, [0 1]);
            xalpha = double(inside);
            xalpha(border) = alpha_border;

            % final smoothing of alpha
            if X.controls.final__alpha__smooth
                xalpha = brick.filt(xalpha,3,'l',[1 2]);
                e = 1e-2;
                xalpha(xalpha<e) = 0;
                xalpha(xalpha>1-e) = 1;
                % -> this extends 'border'
                inside = (xalpha == 1);
                outside = (xalpha == 0);
                border = ~inside & ~outside;
                % repeat previous computation
                input_border = brick.imvect(X.input,border,'vector');
                darkness = 1 - mean(input_border,2);
                truecolor_border = inside_spread(border, :);
                truecolor_border_darkness = 1 - mean(truecolor_border,2);
                alpha_border = brick.clip(darkness ./ truecolor_border_darkness, [0 1]);
            end
            xtruecolor(outside, :) = 1;

            % True color can be improved where alpha is large enough so we
            % can trust the pixel color
            %  we have:    input = alpha * truecolor + (1-alpha)
            %  hence:      truecolor = 1 - (1 - input)/alpha
            truecolor_border2 = brick.clip(1 - brick.div(1-input_border,alpha_border), [0 1]);
            truecolor_border2(alpha_border==0, :) = 1;
            smooth = X.controls.true__color__smoothing;
            if ~ismember(smooth, [0 1])
                smooth = alpha_border .^ atanh(1-smooth);
            end
            xtruecolor(border, :) = brick.mult(smooth,truecolor_border) + brick.mult(1-smooth,truecolor_border2);

            % true color
            xtruecolor = brick.imvect(xtruecolor,[X.nx X.ny],'image');

            % save result, possibly only for sub-part of image
            X.save_result(xalpha, xtruecolor);
            
            % display result
            X.display_final_result()
        end
        function save_result(X, xalpha, xtruecolor)
            if ~isempty(X.dosubimage)
                X.truecolor = brick.imvect(X.truecolor);
                xtruecolor = brick.imvect(xtruecolor);
                switch X.dosubimage
                    case 'inside'
                        mask = X.image_subpart(:);
                    case 'outside'
                        mask = ~X.image_subpart(:);
                end
                X.truecolor(mask,:) = xtruecolor(mask,:);
                X.truecolor = brick.imvect(X.truecolor,[X.nx X.ny]);
                X.alpha(X.image_subpart) = xalpha(X.image_subpart);
            else
                X.truecolor = xtruecolor;
                X.alpha = xalpha;
            end
        end
        function display_intermediary_result(X, flag)
            if nargin>=2 && strcmp(flag, 'original')
                raw = X.original_input;
            else
                raw = X.input;
            end
            set(X.im.raw,'cdata',permute(raw,[2 1 3]))
            set(X.im.masks,'cdata',permute(X.borders,[2 1 3]))
        end
        function display_final_result(X)
            b = brick.add(brick.mult(X.truecolor,X.alpha), X.checker.*(1-X.alpha));
            set(X.im.truecolor,'cdata',permute(X.truecolor,[2 1 3]))
            set(X.im.alpha,'cdata',permute(X.alpha,[2 1 3]))
            set(X.im.result,'cdata',permute(b,[2 1 3]))
        end
        function save(X)
            % file name
            if isempty(X.filename)
                fsave = '*.png';
            else
                fsave = [brick.fileparts(X.filename,'base') ' - transparency.png'];
            end
            fsave = brick.savefile(fsave,'Save image with transparency as');

            % save image
            xtruecolor = X.truecolor;
            if size(xtruecolor,3)==1
                xtruecolor = repmat(xtruecolor,[1 1 3]);
            end
            a = cat(3,xtruecolor,X.alpha);
            brick.saveimg(a,fsave)
        end
        function copy_clipboard(X)
            xtruecolor = X.truecolor;
            if size(xtruecolor,3)==1
                xtruecolor = repmat(xtruecolor,[1 1 3]);
            end
            a = cat(3,xtruecolor,X.alpha);
            imclipboard('copy', permute(a,[2 1 3]))
        end
        function paste_clipboard(X)
            a = permute(imclipboard('paste'),[2 1 3]);
            X.load_image(a)
        end
    end
    
    

end