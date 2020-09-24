% load data
data = load(fullfile(xplor.demo.demo_data_folder, 'intrinsic.mat'));

% brain activity data measured using intrinsic imaging
x = xplr.XData(data.x, ...
    {{'x' 'um' 24}, {'y' 'um' 24}, {'time' 's' 0.2}, {'condition' {'stimulated' 'control'}}, 'repetition'}, ...
    'INTRINSIC IMAGING: WHISKER STIMULATION');

% image display of activity
xplor(x, ...
    'view&ROI', {'x' 'y'}, ...
    'colormap', 'gray')

% time display of activity
xplor(x, ...
    'view&ROI', 'time')

% display of reference anatomic image
xplor(data.img, ...
    'header', {{'x' 'um' 4}, {'y' 'um' 4}}, ...
    'name', 'BRAIN AREA', ...
    'view&ROI', {'x' 'y'}, ...
    'colormap', 'gray')

% display some info
disp('---')
disp('Intrinsic Imaging demo:')
disp('Drag in the displays with the mouse left or right button')
disp('to select what is shown in other displays.')
disp('---')
