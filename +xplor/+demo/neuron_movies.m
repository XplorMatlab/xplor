% Data
s = brick.loadvar(fullfile(xplr.demo.demo_data_folder,'neuron_movies.mat'));

% Add header information
x = xplr.XData(s.data, ...
    {{'x' 'um' s.pixel_size}, {'y' 'um' s.pixel_size}, {'time' 's' 1/s.frame_rate}, 'trial'}, ...
    'NEURONS');

% Spatial display
V1 = xplor(x, ...
    'view&ROI', {'x' 'y'});
V1.hf.Position = [24   389   560   473];
V1.D.color_map.c_map_def = 'gray';
V1.D.navigation.cross_color = 'white';

% Time display
V2 = xplor(x, ...
    'view&ROI', 'time');
V2.hf.Position = [591   370   685   645];

