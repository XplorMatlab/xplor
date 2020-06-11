
%% Download data
disp 'Load coronavirus data from opendata.ecdc.europa.eu'
% url = 'https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide.csv';
url = 'https://opendata.ecdc.europa.eu/covid19/casedistribution/csv';
file = fullfile(fileparts(which('xplor')),'demo','coronavirus.csv');
websave(file,url);
data = readtable(file);

%% Header information

disp 'Process data'

% Dates
dates_row = data(:,1); % first column
dates_row = table2array(dates_row); % convert to datetime class
date_start = min(dates_row);
dates_row_num = days(dates_row - date_start); % convert to numeric, starting from day 0
nday = max(dates_row_num) + 1;
dates = date_start + (0:nday-1);
dates = fn_map(dates,@(t)datestr(t,'dd/mm'));

% Countries
countries = table2array(data(:,7)); % cell of char arrays
[countries, idxcountry2row, idxrow2country] = unique(countries);
countries = strrep(countries,'_',' ');
pop = table2array(data(idxcountry2row,10))';
ncountry = length(countries);

%% Data

% new cases and deaths
newcasesanddeaths = table2array(data(:,[5 6]));
CORONAVIRUS = zeros(nday, ncountry, 2);
for i = 1:size(data,1)
    CORONAVIRUS(dates_row_num(i)+1, idxrow2country(i), :) = newcasesanddeaths(i,:);
end

% normalize by population
CORONAVIRUS = cat(4,CORONAVIRUS,fn_div(CORONAVIRUS,pop/1e6));

% cumulated sums
CORONAVIRUS = cat(5,CORONAVIRUS,cumsum(CORONAVIRUS,1));

datanames = {'cases' 'deaths'; ...
    'total' 'per million people'; ...
    'dayly' 'cumulated'};

% % normalize by population
% CORONAVIRUS = cat(3,CORONAVIRUS,fn_div(coronavirus,pop));
% % cumulated sums
% CORONAVIRUS = cat(3,CORONAVIRUS,cumsum(coronavirus,1));
% % deads/cases
% CORONAVIRUS = cat(3,CORONAVIRUS,CORONAVIRUS(:,:,6)./CORONAVIRUS(:,:,5));
% datanames = {'new cases' 'new deaths' 'new cases/population' 'new deaths/population' ...
%     'cases' 'deaths' 'cases/population' 'deaths/population' ...
%     'deaths/cases'};

sdata = size(CORONAVIRUS); sdata = sdata(3:end);

%% Get world map

folder = fullfile(fileparts(which('xplor')),'demo');
world_map_mat = fullfile(folder,'world_map.mat');
if exist(world_map_mat, 'file')
    disp 'load world map'
    load(world_map_mat)
else
    %% Get world map from internet (requires Mapping Toolbox and Internet connection)

    folder = fullfile(fileparts(which('xplor')),'demo');
    worldmapfile = fullfile(folder,'worldmap.shp');
    if ~exist(worldmapfile,'file')
        %% ()
        disp 'Get world map from thematicmapping.org'
        % see https://thematicmapping.org/downloads/world_borders.php
        url = 'https://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip';
        files = unzip(url, folder);
        for i = 1:length(files)
            f = files{i};
            if strcmp(fn_fileparts(f,'name'), 'Readme.txt')
                delete(f)
            else
                movefile(f, strrep(f, 'TM_WORLD_BORDERS_SIMPL-0.3', 'worldmap'))
            end
        end        
    end
    S = shaperead(worldmapfile);
    nshape = length(S);

    %% Apply the Robinson projection to all country shapes

    % projection structure
    mstruct = defaultm('robinson'); 
    mstruct.origin = [0 0 0]; 
    mstruct = defaultm(mstruct);

    % bounding box
    [left, ~]   = mfwdtran(mstruct,0,-180);
    [right, ~]  = mfwdtran(mstruct,0,180);
    [~, bottom] = mfwdtran(mstruct,-90,0);
    [~, top]    = mfwdtran(mstruct,90,0);
    boundingbox = [left bottom; right top]; % columns = X/Y

    % apply projection
    for i = 1:nshape
        % area-preserving projection
        [x, y] = deal(S(i).X, S(i).Y);                  % longitude/latitude
        [S(i).x, S(i).y] = mfwdtran(mstruct,y,x);       % apply projection
    end

    %% Save
    save(world_map_mat,'S','boundingbox')
end

%% Project map shapes into image

% image size and scaling
nx = 1000; % image width
scale = (nx-2.2) / diff(boundingbox(:,1));
ny = ceil(diff(boundingbox(:,2)) * scale + 1.7);
scale(2,1) = -scale; % Y goes from bottom to top, but image indices go from top to bottom
offset = [1.6; ny-.6] - boundingbox(1,:)'.*scale;     % offset(1) + bounddingbox(2,1) will be equal to [1+0.6 nx-0.6]

% country colors
nshape = length(S);
colors = randomcolors(nshape);

% make map!
map = zeros(nx,ny);
mapcolor = zeros(nx*ny, 3);
for i = 1:nshape
    % area-preserving projection
    [x, y] = deal(S(i).x, S(i).y);        % coordinates in the Robinson projection
    [x, y] = deal(offset(1) + scale(1)*x, offset(2) + scale(2)*y);  % coordinates in image
    mask = fn_poly2mask(x, y, nx, ny);
    np = sum(mask(:));
    map(mask) = i;
    mapcolor(mask, :) = repmat(colors(i,:), [np 1]);
end
mapcolor = reshape(mapcolor,[nx ny 3]);

% % display
% figure(1)
% imagesc(permute(mapcolor,[2 1 3]))
% axis image

%% Link coronavirus dataset countries with map shapes

shapenames = lower({S.NAME});

country2shape = zeros(1,ncountry);

% look for matches
for i = 1:ncountry
    name = lower(countries{i});

    % some tiny differences
    name = fn_switch(name, ...
        'cote divoire',     'cote d''ivoire', ...
        'laos',             'lao', ...
        'timor leste',      'timor', ...
        'north macedonia',  'macedonia', ...
        'south korea',      'korea, republic of', ...
        'united states of america',     'united states', ...
        'vietnam',          'viet nam', ...
        name);
    
    idx = fn_find(name,shapenames,'first');
    if isempty(idx)
        idx = fn_find(@(s)strfind(s,name), shapenames, 'first');
    end
    if ~isempty(idx)
        country2shape(i) = idx;
    else
        disp(['did not find ' lower(countries{i})])
    end
end

% remove countries not found in shapes
unfound = (country2shape == 0);
countries(unfound) = [];
country2shape(unfound) = [];
ncountry = length(countries);
CORONAVIRUS(:,unfound,:,:,:) = [];

%% Make movie!!

disp 'Make movie from countries data'
CORONAVIRUS_MAP = NaN(nday, nx*ny, prod(sdata));
population = NaN(nx,ny);
for i = 1:ncountry
    mask = (map == country2shape(i));
    CORONAVIRUS_MAP(:, mask, :) = repmat(CORONAVIRUS(:, i, :), [1 sum(mask(:)) 1]);
    population(mask) = pop(i);
end
CORONAVIRUS_MAP = reshape(CORONAVIRUS_MAP, [nday nx ny sdata]);

%% Save data

if eval('false')
    % ()
    disp 'Save'
    cd(folder)
    save coronavirus CORONAVIRUS CORONAVIRUS_MAP dates nday countries ncountry datanames pop population
end