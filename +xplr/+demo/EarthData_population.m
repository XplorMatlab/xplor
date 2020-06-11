%% Download data from (needs registration):
% https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11/data-download

%% Load data

% resolution 1440x720 (data at original resolution 8640x4320 can be found
% on the website)
folder = fullfile(demo_data_folder,'EarthData','gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11_totpop_15_min_nc');
file = fullfile(folder,'gpw_v4_population_count_adjusted_rev11_15_min.nc');
ncdisp(file)
longitude = ncread(file,'longitude');
latitude = ncread(file,'latitude');
years = {'2000' '2005' '2010' '2015' '2020'};
data = ncread(file,'UN WPP-Adjusted Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 15 arc-minutes');
population = data(:,:,1:5);
country = data(:,:,11); country(country==max(country(:))) = NaN;

%% XPLOR it
x = xplr.xdata(population, ...
    {{'longitude' '°' diff(longitude(1:2)) longitude(1)}, {'latitude' '°' diff(latitude(1:2)) latitude(1)}, {'year' years}}, ...
    'POPULATION');

%% 
xplor x