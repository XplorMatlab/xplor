% This animation illustrates the idea of navigation inside a
% simple spatio-temporal data, but does not actually use XPLOR!

%% First idea: image engine.jpg

% %% Load engine picture
% file = fullfile(xplor.demo.demo_data_folder, 'engine.jpg');
% engine_color = brick.readimg(file);
% engine_color = double(engine_color) / 255;
% 
% %% Analyze colormap to convert color image back to intensities image
% brick.figure('engine', 'color', 'w')
% brick.show_image(engine_color)
% cmap = fliplr(engine_color([372:375 388:390],16:275,:));
% cmap = squeeze(mean(cmap, 1));
% brick.showcolormap(cmap)
% engine = brick.color2index(engine_color, cmap);
% engine = (engine - 1) / (length(cmap) - 1);
% 
% % crop
% engine = engine(43:320,:);
% 
% % display
% figure(1)
% subplot(111)
% brick.imvalue('image')
% imagesc(engine', [0 1]), axis image
% colormap(cmap)
% 
% %% Decompose image between "base heating" and "over heating"
% thr = .3;
% alpha = .1;
% xbase = min(engine, thr); xbase = xbase + alpha * (engine-xbase);
% mod = engine - xbase;
% 
% % display
% subplot(111)
% imagesc(xbase', [0 1]), axis image
% colormap(cmap)
% 
% %% Time
% [nx, ny] = size(engine);
% T = 100; % 10 instants
% idx = 1:10:T;
% ni = length(idx);
% 
% %% Heat field
% sigmax = 200; % spatial smoothing
% sigmat = 30; % temporal smoothing
% noise_level = 1; % std of noise
% field = randn(nx, ny, T + 4*ceil(sigmat));
% field = brick.filt(field, sigmax, 'l', [1 2]);
% field = brick.filt(field, sigmat, 'l', 3, 'pink');
% field = field(:, :, 2*ceil(sigmat)+(1:T));
% field = max(0, (field - prctile(field, 20)));
% field = field / std(field(:)) * .5;
% 
% % display
% figure(1), subplot(211)
% brick.framedisplay(field(:,:,idx), 'fit', 'singleaxes', 'ncol', ni)
% colormap(cmap)
% 
% subplot(212)
% x = xbase + mod.*field;
% brick.framedisplay(x(:,:,idx), [0 1], 'singleaxes', 'ncol', ni)

% This animation illustrates the idea of navigation inside a
% simple spatio-temporal data, but does not actually use XPLOR!

%% Load engine picture
file = fullfile(xplor.demo.demo_data_folder, 'engine.png');
engine = brick.readimg(file);
engine = double(engine) / 255;
engine = brick.bin(engine, 2);
mask = engine(:,:,4);
engine = mean(engine(:,:,1:3), 3);

brick.figure('engine', 'color', 'w')
brick.show_image(engine)

%% Color map: mapgeog
cmap = colormaps.mapgeog(256);

%% Time
[nx, ny] = size(engine);
T = 100; % 10 instants
idx = 1:10:T;
ni = length(idx);

%% Heat field
sigmax = 200; % spatial smoothing
sigmat = 30; % temporal smoothing
noise_level = 1; % std of noise
field = randn(nx, ny, T + 4*ceil(sigmat));
field = brick.filt(field, sigmax, 'l', [1 2]);
field = brick.filt(field, sigmat, 'l', 3, 'pink');
field = field(:, :, 2*ceil(sigmat)+(1:T));
field = max(0, (field - prctile(field, 20)));
field = field / std(field(:)) * .5;

%% display
figure(1), subplot(211)
brick.framedisplay(field(:,:,idx), 'fit', 'singleaxes', 'ncol', ni)
colormap(cmap)

subplot(212)
x = engine/2 + mask.*field/2;
brick.framedisplay(x(:,:,idx), [0 1], 'singleaxes', 'ncol', ni)

