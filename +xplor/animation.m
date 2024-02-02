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
if ~exist('engine', 'var')
    file = fullfile(xplor.demo.demo_data_folder, 'engine.png');
    engine = brick.readimg(file);
    engine = double(engine) / 255;
    engine = brick.bin(engine, 2);
    mask = engine(:,:,4);
    engine = mean(engine(:,:,1:3), 3);

    brick.figure('engine', 'color', 'w')
    brick.show_image(engine)
end

%% Time
[nx, ny] = size(engine);
r = nx/ny;
T = 100; % 10 instants
idx = 1:20:T;
ni = length(idx);

%% Color map: hot
cmap = hot(256);
amap = max(cmap, [], 2);
cmap2 = (1-amap) * .5 + cmap .* amap; % apply transparency with gray background
amap = amap * .75;

%% Heat field
if ~exist('field', 'var')
    rng(5)
    sigmax = 200; % spatial smoothing
    sigmat = 100; % temporal smoothing
    noise_level = 1; % std of noise
    field = randn(nx, ny, T + 4*ceil(sigmat));
    field = brick.filt(field, sigmax, 'l', [1 2]);
    field = brick.filt(field, sigmat, 'l', 3, 'pink');
    field = field(:, :, 2*ceil(sigmat)+(1:T));
    field = max(0, (field - prctile(field, 20)));
    field = field / std(field(:)) * .4;

    clip = [0 1.5];
    movie = brick.clip(field, clip, cmap);
    alpha = brick.clip(field, clip, amap);
    merge = .75;
    movie = (1-alpha) .* engine + alpha .* movie;
    movie_masked = (1-mask) + mask.*movie;
end

%% display
if eval('false')
    %%
    figure(1), subplot(211)
    brick.framedisplay(field(:,:,idx), 'fit', 'singleaxes', 'ncol', ni)
    colormap(cmap2)

    subplot(212)
    clip = [0 1.5];
    brick.framedisplay(movie_masked(:,:,:,idx), [0 1], 'singleaxes', 'ncol', ni)
end

%% Animation window
W = 800;
H = 400;
R = W/H;
hf = brick.figure('XPLOR animation', [W H], 'color', 'w');
ha = axes('position', [0 0 1 1], 'visible', 'off');
% axis([-.02 1.02 -.02 1.02])
axis([0 1 0 1])
fps = 15;

%% Init movie capture
duration = 2;
nframe = duration * fps;
capture = zeros(W, H, 3, 4*nframe, 'uint8');

%% Part 1: Movie
figure(hf)
cla(ha)

% prepare 5 images
nimg = 5;
im = gobjects(1, nimg);
x1 = .5+.5*[-1; 1]*r/R;
y1 = [1 0];
for k = 1:nimg
    im(k) = surface( ...
        [x1 x1], [y1; y1], zeros(2), ...
        'parent', ha, ...
        'EdgeColor', 'none', 'FaceColor', 'texturemap', ...
        'CDataMapping', 'scaled', 'FaceAlpha', 'texturemap', ...
        'CData', movie(:,:,:,1), 'AlphaData', mask);
end
im = fliplr(im);
set(im(2:end), 'visible', 'off')

step = floor(T / nframe);
nperblock = 2*fps / nimg;
for i = 1:nframe
    tic
    kfr = 1 + (i-1)*step;
    set(im(1), 'cdata', movie(:,:,:,kfr))
    M = getframe(hf);
    capture(:,:,:,i) = permute(M.cdata, [2 1 3]);
    pause(1/fps - toc)
end

%% Part 2: Spread movie

% max-movie
% dampen_max = .7;
% binned_field = brick.bin(field(:,:,1:nimg*nperblock*step), [1 1 nperblock*step], 'max') * dampen_max;
binned_field = brick.bin(field(:,:,1:nimg*nperblock*step), [1 1 nperblock*step]);
binned_movie = brick.clip(binned_field, clip, cmap);
binned_alpha = brick.clip(binned_field, clip, amap);
binned_movie = (1-binned_alpha) .* engine + binned_alpha .* binned_movie;

overlap = .35;
y_image = .5;
w = 1 / (nimg - (nimg-1)*overlap);
dx = w * (1-overlap);
h = R/r*w;
for k = 1:nimg
    x2 = (k-1) * dx + [0; w];
    y2 = y_image + .5*[1 -1]*h;
    offset = 1 + (k-1)*nperblock*step;
    for i = 1:nperblock
        tic
        x = x1 + i/nperblock * (x2-x1);
        y = y1 + i/nperblock * (y2-y1);
        kfr = offset + (i-1)*step;
        set(im(k),...
            'visible', 'on', ...
            'xdata', [x x], 'ydata', [y; y], ...
            'cdata', movie(:,:,:,kfr))
        M = getframe(hf);
        capture(:,:,:,nframe+(k-1)*nperblock+i) = permute(M.cdata, [2 1 3]);
        pause(1/fps - toc)        
    end
    %     i = 1 + (k-1)/(nimg-1)*(nperblock-1);
    %     kfr = round(offset + (i-1)*step);
%     set(im(k), 'cdata', movie(:,:,:,kfr))
    set(im(k), 'cdata', binned_movie(:,:,:,k))
    [x1, y1] = deal(x2, y2);
end

%% Part 3: ROIs and time courses

figure(hf)
delete(findobj('type', 'line'))

ROI = [270 90; 330 150; 330 270]; % steps of 60
radius = 30;
colors = {[0 0 1], [0 1 0], [1 .1 .1]};
graph_y0 = [.67 .42 .18];
% graph_y00 = -.05;
amplitude = .1;

nroi = size(ROI, 1);

nsub = 2; % subsampling
ntgraph = 1 + (nimg*nperblock-1)*step;
ntgraphsub = 1 + (nimg*nperblock-1)*step*nsub;
graph_data = zeros(ntgraphsub, nroi);
graph_xstart = (w-dx)/2;
graph_xstop = 1 - graph_xstart;
graph_y = NaN(ntgraphsub, nroi);
graph_line = gobjects(1, nroi);
alpha3 = zeros(nx, ny);
for j = 1:nroi
    % ROI
    sel = xplr.SelectionND('ellipse2D', {ROI(j,:) radius 1}, [nx ny]);
    alpha3(sel.data_ind) = 1;

    % extract graph data
    xi = brick.roiavg(field, sel.data_ind);
    xi = xi(1:1+(nimg*nperblock-1)*step);
    graph_data(:, j) = interp1(xi, 1:1/nsub:ntgraph(end), 'spline');
    
    % mark ROI
    for k = 1:nimg
        x2 = (k-1) * dx + [0; w];
        y2 = y_image + .5*[1 -1]*h;
        poly = sel.polygon; % image coordinates
        poly = (poly-.5) ./ [nx; ny];  % normalized coordinates
        poly = [x2(1); y2(1)] + [diff(x2); diff(y2)].*poly; % axes coordinates
        brick.drawpoly(poly, 'color', colors{j}, 'linewidth', 2)
    end
    
    % position graph line
    %     x1start = [0; w];
    %     x1stop = (nimg-1) * dx + [0; w];
    %     xroi = (ROI(j,1)-.5)/nx;                % normalized coordinates
    %     xstart = x1start(1) + diff(x1start)*xroi; % axes coordinates
    %     xstop = x1stop(1) + diff(x1stop)*xroi; % axes coordinates
%     y1 = y_image + .5*[1 -1]*h;
%     yroi = (ROI(j,2)-.5)/ny;                  % normalized coordinates
%     graph_y0(j) = y1(1) + diff(y1)*yroi + graph_y00; % axes coordinates
    graph_line(j) = line(linspace(graph_xstart,graph_xstop,ntgraphsub), ...
        graph_y(:, j), ...
        'color', colors{j}, 'linewidth', 2);
end
movie3 = (1-alpha3).*engine + alpha3.*movie;

% Make sure images are at the right position
for k = 1:nimg
    x2 = (k-1) * dx + [0; w];
    y2 = y_image + .5*[1 -1]*h;
    set(im(k), 'xdata', [x2 x2], 'ydata', [y2; y2], ...
        'cdata', binned_movie(:,:,:,k), 'visible', 'on')
end

% play movie in 5 subblocks and draw graphs
for k = 1:nimg
    x2 = (k-1) * dx + [0; w];
    y2 = y_image + .5*[1 -1]*h;
    for i = 1:nperblock
        tic
        x = x1 + i/nperblock * (x2-x1);
        y = y1 + i/nperblock * (y2-y1);
        kfr = 1 + ((k-1)*nperblock+(i-1))*step;
        set(im(k), 'cdata', movie3(:,:,:,kfr))
        for j = 1:nroi
            idx = 1:1+(kfr-1)*nsub;
            graph_y(idx, j) = graph_y0(j) + amplitude * graph_data(idx, j);
            set(graph_line(j), 'ydata', graph_y(:, j))
        end
        M = getframe(hf);
        capture(:,:,:,2*nframe+(k-1)*nperblock+i) = permute(M.cdata, [2 1 3]);
        pause(1/fps - toc)        
    end
%     i = 1 + (k-1)/(nimg-1)*(nperblock-1);
%     kfr = round(1 + ((k-1)*nperblock+(i-1))*step);
%     set(im(k), 'cdata', movie(:,:,:,kfr))
    set(im(k), 'cdata', (1-alpha3).*engine + alpha3.*binned_movie(:,:,:,k))
    [x1, y1] = deal(x2, y2);
end


%% Part 4: Many ROIs arranged as a grid

figure(hf)
delete(findobj('type', 'line'))

% Create ROIs and extract their signals
[xrois, yrois] = ndgrid(30:60:nx, 30:60:ny);
ok = false(size(xrois));
alpha4 = zeros(nx, ny);
graph_data4 = zeros(ntgraphsub, length(xrois), length(yrois));
for i = 1:size(xrois, 1)
    for j = 1:size(xrois, 2)
        sel = xplr.SelectionND('ellipse2D', {[xrois(i,j) yrois(i,j)] radius 1}, [nx ny]);
        ok(i,j) = (mean(mask(sel.data_ind) > .5) > .1); % at least 10% of the ROI on the engine
        if ~ok(i,j), continue, end
        alpha4(sel.data_ind) = 1;
        xij = brick.roiavg(field, sel.data_ind);
        xij = xij(1:1+(nimg*nperblock-1)*step);
        graph_data4(:, i, j) = interp1(xij, 1:1/nsub:ntgraph(end), 'spline');
    end
end
xrois = xrois(ok);
yrois = yrois(ok);
graph_data4 = graph_data4(:, ok);
nroi4 = length(xrois);

% Movie: show signal only inside the ROIs
movie4 = (1-alpha4).*engine + alpha4.*movie;
dim0 = .6; % make image dimmer
dim = dim0;

% Positions for displaying ROI signals
gap = .25; % horizontal gap between signals
ntshow4 = 35; % temporal zooming
yoffset0 = -.1;
yamplitude4 = .13;
centers = [xrois(:)'; yrois(:)']; % image coordinates
centers = (centers - .5) ./ [nx; ny]; % normalized coordinates
centers(2,:) = 1 - centers(2,:);      % axes coordinates (image first row is at the top)
w4 = (1-gap) * (60/nx); % width
xdata = centers(1, :) + linspace(-w4/2, w4/2, ntshow4)';
yoffset_points = yoffset0 + centers(2, :);
yoffset = yoffset_points + zeros(ntshow4, nroi4);
graph_line4 = line(xdata, yoffset, 'color', 'k', 'linewidth', 1.2);
graph_point4 = gobjects(1, nroi4);
for j = 1:nroi4
    graph_point4(j) = line(centers(1,j), yoffset_points(j), 'color', 'k', 'linestyle', 'none', ...
        'marker', '.', 'markersize', 25);
end

% Keep colors of the 3 previous ROIs
for j = 1:nroi
    idx = find(xrois==ROI(j,1) & yrois==ROI(j,2));
    set(graph_line4(idx), 'color', colors{j}, 'linewidth', 2)
    set(graph_point4(idx), 'color', colors{j}, 'linewidth', 2)
end

% Images
set(im(2:end), 'visible', 'off')  % show only first one
x2 = (0:nimg-1) * dx + [0; w];
y2 = y_image + .5*[1 -1]*h;
x4 = [0; 1];
y4 = [1 0];
x1 = .5+.5*[-1; 1]*r/R;
set(im, 'visible', 'on')
for k = 1:nimg
    set(im(k), 'xdata', x2(:, [k k]), 'ydata', [y2; y2])
end

% Play movie
% transition
nonset = 1;  % here use value 1 (not 0) for no transition
noffset = 0; % here use value 0 for no transition
for i = 1:nframe
    tic
    
    % current frame
    kfr = 1 + (i-1)*step;
    
    % move images
    if i <= nonset
        x = x2 + (i/nonset)*(x4 - x2);
        y = y2 + (i/nonset)*(y4 - y2);
        for k = 1:nimg
            set(im(k), 'xdata', x(:, [k k]), 'ydata', [y; y])
        end
        if i==nonset
            set(im(2:4), 'visible', 'off')
        end
    elseif i >= nframe-(noffset-1)
        ratio = r/R + (nframe-i)/noffset * (1-r/R);
        set(im, 'xdata', .5+([x4 x4]-.5)*ratio)
        dim = 1 + (nframe-i)/noffset * (dim0 - 1);
    end
    
    % image data
%     x = x0 + i/nperblock * (x1-x0);
%     y = y0 + i/nperblock * (y1-y0);
%     kfr = 1 + ((k-1)*nperblock+(i-1))*step;
    set(im, 'cdata', (1-dim) + dim*movie4(:,:,:,kfr))


    % graph data
    kfrsub = 1 + (kfr-1)*nsub;
    frsubidx = kfrsub + (-floor(ntshow4/2):floor(ntshow4/2));
%     frsubidx = frsubidx + max(0, 1-frsubidx(1)) + min(0, ntgraphsub-frsubidx(end)); % make sure values between 1 and ntgraphsub
    okfr = (frsubidx >=1 & frsubidx <= ntgraphsub); 
    yi = NaN(ntshow4, nroi4);
    yi(okfr, :) = yoffset(okfr, :) + yamplitude4 * graph_data4(frsubidx(okfr), :);
    for j = 1:nroi4
        set(graph_line4(j), 'ydata', yi(:, j))
        set(graph_point4(j), 'ydata', yoffset_points(j) + yamplitude4 * graph_data4(kfrsub, j))
        if i >= nframe-(noffset-1)
            set(graph_line4(j), 'xdata', .5 + (xdata(:,j)-.5)*ratio)
            set(graph_point4(j), 'xdata', .5 + (centers(1,j)-.5)*ratio)
        end
    end
    M = getframe(hf);
    capture(:,:,:,3*nframe+i) = permute(M.cdata, [2 1 3]);
    pause(1/fps - toc)
end

%% Part 1: Movie [repeated to see the transition]

if eval('false')
    figure(hf)
    cla(ha)
    
    % prepare 5 images
    nimg = 5;
    im = gobjects(1, nimg);
    x1 = .5+.5*[-1; 1]*r/R;
    y1 = [1 0];
    for k = 1:nimg
        im(k) = surface( ...
            [x1 x1], [y1; y1], zeros(2), ...
            'parent', ha, ...
            'EdgeColor', 'none', 'FaceColor', 'texturemap', ...
            'CDataMapping', 'scaled', 'FaceAlpha', 'texturemap', ...
            'CData', movie(:,:,:,1), 'AlphaData', mask);
    end
    im = fliplr(im);
    set(im(2:end), 'visible', 'off')
    
    step = floor(T / nframe);
    nperblock = 2*fps / nimg;
    for i = 1:nframe
        tic
        kfr = 1 + (i-1)*step;
        set(im(1), 'cdata', movie(:,:,:,kfr))
        M = getframe(hf);
        capture(:,:,:,i) = permute(M.cdata, [2 1 3]);
        pause(1/fps - toc)
    end
end

%% Save movie
cd('G:\.shortcut-targets-by-id\17q2UWinYzrqTTrLISieZjykHY16apTm-\Learning Robots - Private\5 - Business\2021-10 Maturaction\Maturaction 2021 - Projet B - Logiciels data-IA\Groupe XPLOR')
brick.savemovie(capture, 'animation XPLOR', 'fps', 15)