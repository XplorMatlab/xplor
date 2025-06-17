%% Load data

files = {'673nmos', '656nmos', '502nmos'};
folder = fullfile(fileparts(which('xplor')),'demo data','Hubble');
brick.mkdir(folder)

% Download from internet
if ~exist(fullfile(folder,[files{1} '.fits']),'file')
    %% (download data from Hubble website
    disp 'Download Eagle Nebula image from internet'
    base_url = 'https://www.spacetelescope.org/static/projects/fits_liberator/datasets/eagle/';
    for i = 1:length(files)
        unzip([base_url files{i} '.zip'], folder);
    end
end

% Read files
disp 'Read images'
clear data
for k = 1:length(files)
    fk = fullfile(folder, [files{k} '.fits']);
    info(k) = fitsinfo(fk);
    data(:,:,k) = fitsread(fk);
end
data0 = data;

%% Normalize colors 

data = brick.mult(data0, shiftdim([1 .15 1.3],-1));
% V.data.chgData(data)

%% Header info

labels = {'x' 'y' 'wavelength'};
unit = '''';
table = info(1).PrimaryData.Keywords;
dRAx = table{strcmp(table(:,1),'CD1_1'),2} * 60;
dDECx = table{strcmp(table(:,1),'CD2_1'),2} * 60;
dx = norm([dRAx dDECx]);

step = 3.15 / size(data,1);
wavelengths = {'673nm' '656nm' '502nm'};
sz = size(data);

header = xplr.Header( ...
    {'x'  unit sz(1) dx -sz(1)/2*dx}, ...
    {'y' unit sz(2) dx -sz(2)/2*dx}, ...
    {'wavelength' wavelengths} ...
    );
Eagle_Nebula = xplr.XData(data, header, 'Eagle Nebula');

%% Display

V = xplor(Eagle_Nebula);
V.D.clipping.auto_clip_mode = 'prc.1';
V.D.clipping.adjust_to_view = true;

%% Display some info
disp('---')
disp('Astronomy demo ("Pillars of the creation" in the Eagle Nebula,')
disp('photographed by Hubble space telescope):')
disp('Zoom in the image with mouse left button and by dragging sliders on the sides.')
disp('Change the clipping range using the three small buttons top-right')
disp('(click + or -, click and drag on the black and white button).')
disp('Find more options for the clipping range in the ''Clipping'' menu.')
disp('---')