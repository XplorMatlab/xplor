
url = 'https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide.csv';
data = webread(url);

%% Header information

% Dates
dates = data(:,1); % first column
dates = table2array(dates); % convert to datetime class
dates = days(dates - min(dates)); % convert to numeric, starting from day 0
nday = max(dates) + 1;

% Countries
countries = table2array(data(:,7)); % cell of char arrays
[countries, idxcountry2row, idxrow2country] = unique(countries);
population = table2array(data(idxcountry2row,9));
ncountry = length(countries);

%% Data

cases = table2array(data(:,5));
deaths = table2array(data(:,6));
pop = table2array(data(:,9));
cases_prc_pop = cases ./ pop * 100;
deaths_prc_pop = deaths ./ pop * 100;
deaths_prc_cases = deaths ./ cases * 100;

data_per_row = [cases deaths cases_prc_pop deaths_prc_pop deaths_prc_cases];
datanames = {'cases' 'deaths' 'cases/population' 'deaths/population' 'deaths/cases'};
ndata = length(datanames);

% make a 3D array out of the Excel table
x = zeros(nday, ncountry, ndata);
for i = 1:size(data_per_row,1)
    
    x(dates(i)+1, idxrow2country(i), :) = data_per_row(i,:);
    
end

%% Test

data = rand(72,144);
rv = [0.4 90 0]; 
worldmap('world');
geoshow(data, rv, 'displaytype', 'texturemap');
C = load('coast');
plotm(C.lat, C.long, 'k');


