% Geographical data from the Nasa Earth Observations archive at https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/
%---
% Note: executing this script will not download nor read any data, except
% for the Readme.txt. It will list the data already downloaded. To download
% more data or read downloaded data, execute from inside the "if
% eval('false')" block of the relevant section.

%% Base file and url folders

base_folder = fullfile(fileparts(which('xplor')),'demo','nasa_neo');
fn_mkdir(base_folder);
base_url = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/';

%% Download Readme

file_readme = fullfile(base_folder,'README.txt');
if ~exist(file_readme,'file')
    url_readme = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/README.txt';
    websave(file_readme,url_readme)
end

%% Scan Readme to obtain image types and details

% read Readme.txt
file_readme = fullfile(base_folder,'README.txt');
readme = fn_readtext(file_readme);

% process: separate lines, skip headers
% readme = strsplit(readme,'\n');
nheader = fn_find(@(str)strfind(str,'-------------'),readme,'first');
readme(1:nheader) = [];

% retrieve types and descriptions
nalltype = length(readme);
alltypes = cell(1,nalltype);
alltypes_descriptions = cell(1,nalltype);
for i = 1:nalltype
    [alltypes{i}, alltypes_descriptions{i}] = fn_regexptokens(readme{i},'^([^ ]*) +(.*[^ ]) *$');
end

%% Manual list of datasets with 1 image per month

monthly_types = { ...
    'AMSRE_SSTAn_M'    ... Sea Surface Temperature Anomaly 2002-2011 (1 month - Aqua/AMSR-E)
    'AQUARIUS_SSS_M'   ... Sea Surface Salinity 2011-2015 (1 month)
    'AURA_NO2_M' ...       Nitrogen Dioxide (1 month)   
    'AURA_OZONE_M'     ... Ozone (1 month)
    'AURA_UVI_CLIM_M'  ... UV Index
    'AVHRR_CLIM_M'     ... Average Sea Surface Temperature 1985-1997 (1 month - AVHRR)
    'AVHRR_SST_M'      ... Sea Surface Temperature 1981-2006 (1 month - AVHRR)
    'CERES_INSOL_M'    ... Solar Insolation (1 month)
    'CERES_LWFLUX_M'   ... Outgoing Longwave Radiation (1 month)
    'CERES_NETFLUX_M'  ... Net Radiation (1 month)
    'CERES_SWFLUX_M'   ... Reflected Shortwave Radiation (1 month)
    'GISS_TA_M'        ... Global Temperature Anomaly (1 month)
    'GRACE_LWE_M'      ... Water Equivalent Anomaly 2002-2017 (1 month)
    'MOD10C1_M_SNOW'   ... Snow Cover (1 month - Terra/MODIS)
    'MOD14A1_M_FIRE'   ... Active Fires (1 month - Terra/MODIS)
    'MOD15A2_M_LAI'    ... Leaf Area Index (1 month - Terra/MODIS)
    'MOD17A2_M_PSN'    ... Net Primary Productivity (1 month - Terra/MODIS)
    'MODAL2_M_AER_OD'  ... Aerosol Optical Thickness (1 month - Terra/MODIS)
    'MODAL2_M_AER_RA'  ... Aerosol Particle Radius (1 month - Terra/MODIS, 2005-16)
    'MODAL2_M_CLD_FR'  ... Cloud Fraction (1 month - Terra/MODIS)
    'MODAL2_M_CLD_OT'  ... Cloud Optical Thickness (1 month - Terra/MODIS)
    'MODAL2_M_CLD_RD'  ... Cloud Particle Radius (1 month - Terra/MODIS)
    'MODAL2_M_CLD_WP'  ... Cloud Water Content (1 month - Terra/MODIS)
    'MODAL2_M_SKY_WV'  ... Water Vapor (1 month - Terra/MODIS)
    'MOD_LSTAD_M'      ... Land Surface Temperature Anomaly [Day] (1 month)
    'MOD_LSTAN_M'      ... Land Surface Temperature Anomaly [Night] (1 month)
    'MOD_LSTD_CLIM_M'  ... Average Land Surface Temperature [Day] (1 month)
    'MOD_LSTD_M'       ... Land Surface Temperature [Day] (1 month - Terra/MODIS)
    'MOD_LSTN_CLIM_M'  ... Average Land Surface Temperature [Night] (1 month)
    'MOD_LSTN_M'       ... Land Surface Temperature [Night] (1 month - Terra/MODIS)
    'MOD_NDVI_M'       ... Vegetation Index (1 month - Terra/MODIS)
    'MOP_CO_M'         ... Carbon Monoxide (1 month - Terra/MOPITT)
    'MWOI_SST_M'       ... Sea Surface Temperature 1998+ (1 month MWOI)
    'MY1DMM_CHLORA'    ... Chlorophyll Concentration (1 month - Aqua/MODIS)
    'MYD28M'           ... Sea Surface Temperature (1 month - Aqua/MODIS)
    'MYDAL2_M_AER_OD'  ... Aerosol Optical Thickness (1 month - Aqua/MODIS)
    'MYDAL2_M_AER_RA'  ... Aerosol Particle Radius (1 month - Aqua/MODIS, 2002-16)
    'MYDAL2_M_CLD_FR'  ... Cloud Fraction (1 month - Aqua/MODIS)
    'MYDAL2_M_CLD_OT'  ... Cloud Optical Thickness (1 month - Aqua/MODIS)
    'MYDAL2_M_CLD_RD'  ... Cloud Particle Radius (1 month - Aqua/MODIS)
    'MYDAL2_M_CLD_WP'  ... Cloud Water Content (1 month - Aqua/MODIS)
    'MYDAL2_M_SKY_WV'  ... Water Vapor (1 month - Aqua/MODIS)
    'SWE_M'            ... Snow Water Equivalent (1 month - Passive Microwave, with optical)
    'TRMM_3B43M'       ... Rainfall (1 month - TRMM)
    };
    

%% Download full datasets!

if eval('false')
    %% (go inside this block to execute)
    
    % select datasets to download
    dowload_types = monthly_types;

    % download loop
    for type = dowload_types
        type_str = type{1}; % Get type as a string
        subfolder = fullfile(base_folder,type_str);
        fn_mkdir(subfolder) % Create folder in current Matlab path

        % list from web folder all files belonging to this dataset
        listing = webread([base_url type_str]);
        filenames = regexp(listing,[type_str '[a-zA-Z_\-0-9]*\.FLOAT\.TIFF'],'match');
        if isempty(filenames)
            disp(['no file found for data set ' type_str])
        end

        % download loop
        nfile = length(filenames);
        fn_progress(type_str, nfile)
        for k = 1:nfile
            fn_progress(k)
            filename = filenames{k};
            url = [base_url type_str '/' filename];
            file = fullfile(subfolder, filename);
            if exist(file,'file'), ok = true; continue, end % file already downloaded
            try
                saved = websave(file, url); % Save as files
            catch webread_exception
                error(['url ' url ' not found'])
            end
        end
    end
    
end

%% Download specific files 

if eval('true')
    %% (go inside this block to execute)
    
        % select specific types to download and years range
    types = {'MOD_LSTAD_M' 'MOD_LSTAN_M' 'MOD_LSTD_M' 'MOD_LSTN_M'};
    years_range = 2017:2019;

    % download loop
    for type = types
        type_str = type{1}; % Get type as a string
        subfolder = fullfile(base_folder,type_str);
        fn_mkdir(subfolder) % Create folder in current Matlab path
        n_dates = 1;
        for year = years_range
            for month = 1:12
                file = [type_str '_' num2str(year) '-' num2str(month,'%02d') '.FLOAT.TIFF']; % File formatted as MOP_CO_M_2013-10.TIFF
                url = [base_url type_str '/' file];
                file = fullfile(subfolder, file);
                if exist(file,'file'), continue, end % file already downloaded
                try
                    %data(:,:,n_dates) = webread(strcat(url, file)); % Read
                    saved = websave(file, url); % Save as files
                    disp(['saved ' file])
                catch webread_exception
                    if strcmp(webread_exception.identifier, 'MATLAB:webservices:HTTP404StatusCodeError')
                        disp(['url ' url ' not found, skipping...'])
                        continue;
                    end
                end
                n_dates = n_dates + 1;
            end
        end
    end

end

%% Check years available for each type

% get year range by scanning which first and last file of each folder
alltypes_years_range = struct;
alltypes_empty = false(1,nalltype);
for ktype = 1:nalltype
    type_str = alltypes{ktype}; % Get type as a string
    subfolder = fullfile(base_folder,type_str);
    d = dir(fullfile(subfolder,'*.TIFF'));
    if isempty(d), alltypes_empty(ktype) = true; continue, end
    first_year = str2double(fn_regexptokens(d(1).name,[type_str '_(\d{4})']));
    last_year = str2double(fn_regexptokens(d(end).name,[type_str '_(\d{4})']));
    alltypes_years_range.(type_str) = [first_year last_year];
end
disp(alltypes_years_range)

%% Read data

if eval('false')
    %% (temperature data set, load and save in 'demo data' folder)
    types = {'MOD_LSTD_M'}; % land temperature day
    [data, desc, years_range] = load_data(types);
    cd(fullfile(fileparts(which('xplor')),'demo data'))

    % reorganize data by months
    data = data(:,:,:);
    
    % months    
    nmonth = 12;
    nyear = length(years_range);
    months = cell(nmonth,nyear);
    for i = 1:nmonth
        for j = 1:nyear
            months{i,j} = sprintf('%.2i-%.4i',i,years_range(j));
        end
    end
    months = row(months);
    
    % subselect months
    ok_month = ~isnan(data(265,195,:));
    data = data(:,:,ok_month);
    months = months(:,ok_month);
    
    % compressed version of data: notice that all temperatures are comprised
    % between -25 and 45 degrees
    temperatures = [-25:.5:-19.5 -19:.25:38 38.5:.5:45];
    x = zeros(size(data),'uint8');
    idx = (data<-19);
    x(idx) = round((data(idx)+25)*2) + 1;
    idx = (data>=-19 & data<=38);
    x(idx) = round((data(idx)+19)*4) + 13;
    idx = (data>38);
    x(idx) = round((data(idx)-38)*2) + 241;
    temperatures = [NaN temperatures];
    temperatures_indices = x;
    readme = { ...
        'Earth land temperatures data 2000-2020', ...
        'Nasa Earth Observations dataset (https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/MOD_LSTN_M/', ...
        'See xplor/demo/nasa_neo.m for more general download code.', ...
        '', ...
        'To explore this data:', ...
        '    data = temperatures(1+temperatures_indices);', ...
        '    xplor data', ...
        'Then set dimension labels ''longitude'' ''latitude'', (both with ''°'' unit) and ''months'',', ...
        'and type in scale/values: ''longitude'', ''latitude'', ''months''.'}';
    
    % latitude and longitude
    lat_long_unit = '°';
    longitude = -179.75:.5:179.75;
    latitude = 89.75:-.5:-89.75;
    
    
    fn_savevar('NEO Earth Temperature.mat', ...
        temperatures,temperatures_indices,readme,lat_long_unit,latitude,longitude,months)
    
    %% (custom data set)
    %     types = alltypes(~alltypes_empty);
    %     types = {'MOD10C1_M_SNOW'}; % snow
    [data, desc] = load_data(types);
    
end

function [data, desc, years_range] = load_data(types, xbin, years_range_restrict, months_range)
% function [data, desc, years_range] = load_data(types,[, years_range_restrict[, months_range]])

% get variables from the base space
[base_folder,alltypes,alltypes_descriptions,alltypes_years_range] = dealc( ...
    evalin('base','{base_folder,alltypes,alltypes_descriptions,alltypes_years_range}'));

% optional inputs
if ~exist('xbin','var'), xbin = 5; end
if ~exist('years_range_restrict','var'), years_range_restrict = []; end
if ~exist('months_range','var'), months_range = []; end

% types descriptions
ntype = length(types);
desc = fn_map(@(type)alltypes_descriptions{strcmp(type,alltypes)},types);

% years range
% (code below selects all available years for the selected types)
years_range = fn_map(@(type)alltypes_years_range.(type),types','array');
years_range = [min(years_range(:,1)) max(years_range(:,2))];
% (subselect years range)
if ~isempty(years_range_restrict)
    years_range(1) = max(years_range(1),years_range_restrict(1));
    years_range(2) = min(years_range(2),years_range_restrict(2));
end
fprintf('years range: from %i to %i\n',years_range)
years_range = years_range(1):years_range(2);
nyear = length(years_range);

% subselect months range or keep the default 1:12
if isempty(months_range)
    months_range = 1:12;
end
nmonth = length(months_range);

% spatial binning
nx = 3600/xbin;
ny = 1800/xbin;

% array
data = NaN([nx ny nmonth nyear ntype]);

kfile = 0;
for ktype = 1:ntype
    type_str = types{ktype}; % Get type as a string
    fn_progress(['Load ' desc{ktype}],ntype*nmonth*nyear)
    subfolder = fullfile(base_folder,type_str);
    for kyear = 1:nyear
        year = years_range(kyear);
        for kmonth = 1:nmonth
            month = months_range(kmonth);
            kfile = kfile + 1;
            fn_progress(kfile)
            file = [type_str '_' num2str(year) '-' num2str(month,'%02d') '.FLOAT.TIFF']; % File formatted as MOP_CO_M_2013-10.TIFF
            file = fullfile(subfolder,file);
            if exist(file,'file')
                x = fn_readimg(file);
                % undefined data: replace 99999 by NaN
                x(x==99999) = NaN;
                % bin
                x = fn_bin(x,[-nx -ny]);
                % put in array
                data(:,:,kmonth,kyear,ktype) = x;
            else
                %fprintf('missing data for type %s, month %.2i/%i\n',type_str,month,year);
            end
        end
    end
end

end


