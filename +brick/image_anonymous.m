classdef image_anonymous < brick.interface
    %IMAGE_ANONYMOUS Blur parts of image
    %---
    % function white2alpha(image)
    % function white2alpha([filename|'clipboard'])
    %---
    % Convert image white background to transparency
    % 
    % See also brick.show_image, brick.combine_images, brick.montage,
    % brick.image_anonymous
    properties
        X
        orig_file
        a % original image
        b % modified image
        bb % all steps!
        im
    end
    properties (SetObservable)
        selection_mode = 'rect';
    end
    
    % Init
    methods
        function I = image_anonymous
            hf = figure(923);
            I = I@brick.interface(hf, 'Image Anonymous');
            I.init_grob()
            I.interface_end()
            I.init_parameters()
        end
        function init_grob(I)
            g = struct;
            g.ha = axes;
            g.buttons(1) = uicontrol('string','Get image from clipboard','callback',@(u,e)I.loadImage('clipboard'));
            g.buttons(2) = uicontrol('string','Copy image to clipboard','callback',@(u,e)I.saveImage('clipboard'));
            g.buttons(3) = uicontrol('string','Reset','callback',@(u,e)I.reset());
            g.buttons(4) = uicontrol('string','Undo','callback',@(u,e)I.undo());
            g.buttons(5) = uicontrol('string','Get image from file','callback',@(u,e)I.loadImage());
            g.buttons(6) = uicontrol('string','Save image to file','callback',@(u,e)I.saveImage('auto'));
            g.buttons(7) = uicontrol('string','Bin image (x2)','callback',@(u,e)I.bin(2));
            g.buttons(8) = uicontrol('string','Bin image (x3)','callback',@(u,e)I.bin(3));
            g.X = uipanel();
            I.grob = g;
        end
        function init_parameters(I)
            s = struct('binning', {4 'stepper 1 2 Inf'});
            I.X = brick.control(s, 'in', I.grob.X);
        end
        function init_menus(I)
            init_menus@brick.interface(I)
            m = I.menus.interface;
            brick.propcontrol(I,'selection_mode',{'menuval' 'rect' 'ellipse' 'depends on mouse button'}, ...
                'parent',m,'label','selection mode','separator','on');
            uimenu(I.hf, 'label', 'Help', 'callback', @(u,e)I.show_help, 'separator', 'on')
        end
    end
    
    % Action
    methods
        function loadImage(I, fname)
            if nargin<2
                fname = brick.getfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'}, 'Select Image');
                if ~fname, return, end
            end
            if strcmp(fname, 'clipboard')
                img = imclipboard('paste');
                if isempty(I.im)
                    errordlg 'No image in clipboard'
                    return
                end
                I.orig_file = [];
            else
                img = imread(fname);
                I.orig_file = fname;
            end
            img = brick.float(img, 'image');
            [I.a, I.b] = deal(img);
            I.bb = {I.a};
            I.im = imagesc(I.b,'parent',I.grob.ha);
            axis(I.grob.ha,'image')
            set(I.grob.ha,'xtick',[],'ytick',[])
            set(I.im, 'buttondownfcn', @(u,e)I.anonymous())
            I.display_image()
        end
        function display_image(I)
            set(I.im, 'cdata', I.b)
        end
        function saveImage(I, fname)
            if nargin<2
                fname = brick.savefile({'*.jpg;*.tif;*.png;*.gif','All Image Files'}, 'Select file for saving image');
            elseif strcmp(fname,'auto') && ~isempty(I.orig_file)
                [f, e] = brick.fileparts(I.orig_file,'noext','ext');
                fname = [f '-blur' e];
            end
            if strcmp(fname,'clipboard')
                imclipboard('copy', I.b)
            else
                imwrite(I.b, fname)
            end
        end
        function anonymous(I)
            bin = I.X.binning;
            mode = I.selection_mode;
            if strcmp(mode, 'depends on mouse button')
                if strcmp(get(I.hf, 'selectiontype'), 'normal')
                    mode = 'rect';
                else
                    mode = 'ellipse';
                end
            end
            switch mode
                case 'rect'
                    rect = brick.mouse(I.grob.ha,'rectax-');
                    xstart = max(0,floor(rect(1)/bin)*bin);
                    xend = min(size(I.b,2)-1,ceil(rect(2)/bin)*bin);
                    ystart = max(0,floor(rect(3)/bin)*bin);
                    yend = min(size(I.b,1)-1,ceil(rect(4)/bin)*bin);
                    jj = xstart+1:xend;
                    ii = ystart+1:yend;
                    I.b(ii,jj,:) = brick.bin(I.b(ii,jj,:),bin,@brick.nmean,'same');
                case 'ellipse'
                    ellipse = brick.mouse(I.grob.ha,'ellipse-');
                    ellipse{1} = 1 + (ellipse{1}-1) / bin; % ellipse in the binned image
                    ellipse{2} = ellipse{2} / bin;
                    [ni, nj, nc] = size(I.b);
                    [jj0, ii0] = ellipseInterior(ellipse, ceil([nj ni]/bin));
                    for k = 1:length(jj0)
                        jj = 1+(jj0(k)-1)*bin:min(nj,jj0(k)*bin);
                        ii = 1+(ii0(k)-1)*bin:min(ni,ii0(k)*bin);
                        for c = 1:nc
                            I.b(ii,jj,c) = mean(brick.row(I.b(ii,jj,c)));
                        end
                    end
            end
            I.bb{end+1} = I.b;
            I.display_image()
        end
        function bin(I, factor)
            I.b = brick.bin(I.b, factor);
            I.bb{end+1} = I.b;
            I.display_image()
        end
        function reset(I)
            I.b = I.a;
            I.bb = {I.a};
            I.display_image()
        end
        function undo(I)
            if isscalar(I.bb), return, end
            I.bb(end) = [];
            I.b = I.bb{end};
            I.display_image()
        end
    end
    
    % Help
    methods
        function show_help(I)
            msgbox('hello')
        end
    end
    
end


%---
function [ii jj] = ellipseInterior(poly, sz)

% geometry of the ellipse and quadratic form
c = poly{1};
u = poly{2};
r = norm(u);
u = u/r;
v = [-u(2); u(1)];
U = [u v];
e = poly{3};
dmax = r*max(1,e); % upper bound on distance to center
S = diag([r e*r].^-1);
A = U*S*U';

% mask of the ellipse
imin = max(1,floor(c(1)-dmax));
imax = min(sz(1),ceil(c(1)+dmax));
jmin = max(1,floor(c(2)-dmax));
jmax = min(sz(2),ceil(c(2)+dmax));
[ii jj] = ndgrid(imin:imax,jmin:jmax);
X = [brick.row(ii)-c(1); brick.row(jj)-c(2)];
UtX2 = (U'*X).^2;
XtAX = r^-2 * (UtX2(1,:) + e^-2*UtX2(2,:)); % row vector, = Xk'*A*Xk for each column Xk of X
mask = (XtAX <= 1);
ii = ii(mask);
jj = jj(mask);

end
