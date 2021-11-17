
%% Read data after attempting to update it from the internet
% This data stops 14/12/2020!!! Better read now data from "Our world in
% data" at https://github.com/owid/covid-19-data/tree/master/public/data
% Either https://covid.ourworldindata.org/data/owid-covid-data.csv
% or https://github.com/owid/covid-19-data/blob/master/public/data/ecdc/locations.csv
% + https://github.com/owid/covid-19-data/blob/master/public/data/ecdc/full_data.csv

try
    disp 'Load coronavirus data from opendata.ecdc.europa.eu'
    % url = 'https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide.csv';
    url = 'https://opendata.ecdc.europa.eu/covid19/casedistribution/csv';
    file = fullfile(xplor.demo.demo_data_folder, 'coronavirus.csv');
    websave(file, url);
catch
    disp '(loading latest data from Internet failed, reading data from local file)'
end
data = readtable(file, 'PreserveVariableNames', true);

%% Header information

disp 'Process data'

% Dates
dates_column = data.dateRep;
date_start = min(dates_column);
dates_num = days(dates_column - date_start); % convert to numeric, starting from day 0
nday = max(dates_num) + 1;
dates = date_start + (0:nday-1);

% Countries
countries = data.countriesAndTerritories; % cell of char arrays
[countries, idxcountry2row, idxrow2country] = unique(countries);
countries = strrep(countries,'_',' ');
population = data.popData2019(idxcountry2row)'; % population in Mhab
ncountry = length(countries);

%% Data

% new cases and deaths
newcasesanddeaths = [data.cases data.deaths];
CORONAVIRUS = zeros(nday, ncountry, 2);
for i = 1:size(data,1)
    CORONAVIRUS(dates_num(i)+1, idxrow2country(i), :) = newcasesanddeaths(i,:);
end

% normalize by population
CORONAVIRUS = cat(4, CORONAVIRUS, brick.div(CORONAVIRUS, population/1e6));

% cumulated sums
CORONAVIRUS = cat(5, CORONAVIRUS, cumsum(CORONAVIRUS, 1));

datanames = {'cases' 'deaths'; ...
    'total' 'per million people'; ...
    'dayly' 'cumulated'};

% % normalize by population
% CORONAVIRUS = cat(3,CORONAVIRUS,brick.div(coronavirus,pop));
% % cumulated sums
% CORONAVIRUS = cat(3,CORONAVIRUS,cumsum(coronavirus,1));
% % deads/cases
% CORONAVIRUS = cat(3,CORONAVIRUS,CORONAVIRUS(:,:,6)./CORONAVIRUS(:,:,5));
% datanames = {'new cases' 'new deaths' 'new cases/population' 'new deaths/population' ...
%     'cases' 'deaths' 'cases/population' 'deaths/population' ...
%     'deaths/cases'};

sdata = size(CORONAVIRUS); sdata = sdata(3:end);

%% Get world map

folder = fullfile(fileparts(which('xplor')),'demo data');
world_map_mat = fullfile(folder,'world_map.mat');
if exist(world_map_mat, 'file')
    disp 'load world map'
    load(world_map_mat)
    nshape = length(S);
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
            if strcmp(brick.fileparts(f,'name'), 'Readme.txt')
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

%% Convert map shapes to SelectionND objects

shapes = xplr.SelectionND(nshape);
for i = 1:nshape
    shapes(i) = xplr.SelectionND('poly2D',[S(i).x; S(i).y]);
end

%% Link coronavirus dataset countries with map shapes

shapenames = lower({S.NAME});

country2shape = zeros(1,ncountry);

% look for matches
for i = 1:ncountry
    name = lower(countries{i});

    % some tiny differences
    name = brick.switch_case(name, ...
        'cote divoire',     'cote d''ivoire', ...
        'laos',             'lao', ...
        'timor leste',      'timor', ...
        'north macedonia',  'macedonia', ...
        'south korea',      'korea, republic of', ...
        'united states of america',     'united states', ...
        'vietnam',          'viet nam', ...
        name);
    
    idx = brick.find(name,shapenames,'first');
    if isempty(idx)
        idx = brick.find(@(s)strfind(s,name), shapenames, 'first');
    end
    if ~isempty(idx)
        country2shape(i) = idx;
    else
        % disp(['did not find ' lower(countries{i})])
    end
end

% remove countries not found in shapes
unfound = (country2shape == 0);
countries(unfound) = [];
country2shape(unfound) = [];
ncountry = length(countries);
CORONAVIRUS(:,unfound,:,:,:) = [];
country_shapes = shapes(country2shape);

%% Smooth the data
CORONAVIRUS_RAW = CORONAVIRUS;
CORONAVIRUS = brick.filt(CORONAVIRUS, 15, 'lk', 1);


%% Make movie!!

if eval('false')
    %% ()
    disp 'Make movie from countries data'
    CORONAVIRUS_MAP = NaN(nday, nx*ny, prod(sdata));
    population = NaN(nx,ny);
    for i = 1:ncountry
        mask = (map == country2shape(i));
        CORONAVIRUS_MAP(:, mask, :) = repmat(CORONAVIRUS(:, i, :), [1 sum(mask(:)) 1]);
        population(mask) = population(i);
    end
    CORONAVIRUS_MAP = reshape(CORONAVIRUS_MAP, [nday nx ny sdata]);
end

%% Display data


headers = xplr.Header( ...
    {'day' nday days(1) dates(1)}, ...
    {'country' {'Name', 'ROI2D'} [countries, num2cell(country_shapes)']}, ...
    {'data', {'Name', 'ViewColor'}, [{'cases'; 'deaths'}, {[1 .2 0]; [0 0 0]}]}, ...
    {'norm.' {'total' 'per million people'}}, ...
    {'type' {'daily' 'cumulated'}});
data = xplr.XData(CORONAVIRUS, headers);

V1 = xplor(data, ...
    'view', {'day' 'data' 'type'}, ...
    'display mode', 'time courses');
V1.D.clipping.auto_clip_mode_no_center = 'prc1';
V1.D.clipping.adjust_to_view = true;
V1.D.set_dim_location({'day' 'data' 'type'}, {'x' 'y' 'merged_data'})
V1.D.clipping.set_independent_dim({'data', 'type'})

V2 = xplor(data, ...
    'view', 'country', ...
    'colormap', 'white_red', ...
    'display mode', 'image');
V2.D.clipping.auto_clip_mode_no_center = 'prc.1';
V2.D.clipping.adjust_to_view = false;
V2.D.graph.show_grid_labels = false;
V2.D.graph.show_separation = true;
V2.D.set_dim_location('country', 'xy')

%% Display some info
disp('---')
disp('Coronavirus demo:')
disp('Select in the lists what should be displayed in the main graph.')
disp('You can select multiple values at once, and double-click values to keep them selected.')
disp('Right-click in the list to get a menu with more options.')
disp('The graphs can be re-organized by grabbing and dragging the labels around.')
disp('You can also switch between image and time courses display with the top-right control.')
disp('---')
