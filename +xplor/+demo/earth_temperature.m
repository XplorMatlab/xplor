% Earth temperature data from the Nasa Earth Observations
% See nasa_neo file for download details

% Data (compressed)
s = brick.loadvar(fullfile(xplor.demo.demo_data_folder,'NEO Earth Temperature.mat'));

% Uncompress data
data = s.temperatures(1+s.temperatures_indices);

% Header information
x = xplr.XData(data, ...
    {{'longitude' '°' diff(s.longitude(1:2)) s.longitude(1)}, ...
    {'latitude' '°' diff(s.latitude(1:2)) s.latitude(1)}, ...
    {'months' s.months}}, ...
    'TEMPERATURE');

% Spatial display
V1 = xplor(x, 'view&ROI', {'longitude' 'latitude'});
V1.hf.Position = [11   316   829   364];

% Time display
V2 = xplor(x, 'view&ROI', 'months');
V2.hf.Position = [987   263   712   640];
