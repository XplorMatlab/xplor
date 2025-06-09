%% Load xplor logo
img = brick.readimg(fullfile(fileparts(which('xplor')),'demo data','XPLOR logo.png'));
img = single(img)/255;
logo_rgb = num2cell(img, [1 2]);

%% Add a 3rd dimension with rotated versions of the logo
% coordinates in original image: place white sides outside
white = all(img == 1, 3);
[white_x, white_y] = deal(all(white, 2), all(white, 1));
i_start = find(~white_x, 1, 'first');
i_end = find(~white_x, 1, 'last');
j_start = find(~white_y, 1, 'first');
j_end = find(~white_y, 1, 'last');
if i_end - i_start > j_end - j_start
    [j_start, j_end] = brick.dealc((j_start+j_end)/2 + [-1 1] * (i_end-i_start)/2);
else
    [i_start, i_end] = brick.dealc((i_start+i_end)/2 + [-1 1] * (j_end-j_start)/2);
end
[xx1, yy1] = ndgrid( ...
    interp1([i_start, i_end], [-1, 1], 1:size(img,1), 'linear', 'extrap'), ...
    interp1([j_start, j_end], [-1, 1], 1:size(img,2), 'linear', 'extrap'));

% downsample and rotate
n2 = 200;
% (coordinates in new image)
[xx2, yy2] = ndgrid(linspace(-1, 1, n2), linspace(-1, 1, n2));
xy2 = [brick.row(xx2); brick.row(yy2)];
theta = brick.column(linspace(0, 2*pi, n2));
% (convert to coordinates in the original image where to interpolate)
R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
xy2 = reshape(R * xy2, [n2, 2, n2, n2]);
xx2 = permute(xy2(:, 1, :, :), [3, 4, 2, 1]);
yy2 = permute(xy2(:, 2, :, :), [3, 4, 2, 1]);
% (perform interpolation)
dat = zeros(n2, n2, 3, n2);
for k = 1:3
    dat(:, :, k, :) = interpn(xx1(:, 1), yy1(1, :), logo_rgb{k}, xx2, yy2, 'linear');
end
dat(isnan(dat)) = 1;

%% display
xplor.close_all_windows

header = {{'x', 'px'}, {'y', 'px'}, {'color', {'r', 'g', 'b'}}, {'z', 'px'}};
name = 'LOGO';
data = xplr.XData(dat, header, name);

% by default, all dimensions will be shown, resulting in a grid display of
% mini-images
V1 = xplor(data);
V1.D.set_bin('z', 5, 'private')

% by specifying 'view' option, only dimensions 'x', 'y' and 'color' will
% be displayed (i.e. dimension 'z' will be filtered); these dimensions will
% be organized automatically (x/y dimensions will be shown as image,
% 'color' dimension will be used as a color channel)
V2 = xplor(data, 'view', {'x', 'y', 'color'});

% by specifying 'view' option as a cell array of cell arrays, we also
% specify how the displayed dimensions will be organized: 'z' dimension
% will go on x-axis, 'y' on y-axis, and 'color' on "merged_data", i.e.
% because it has 3 values it will be used as a color channel
V3 = xplor(data, 'view', {{'z'}, {'y'}, {}, {'color'}});






