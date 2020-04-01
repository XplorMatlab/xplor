
%% Download Readme

url_readme = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/README.txt';
folder = fullfile(fileparts(which('xplor')),'demo','nasa_neo');
file_readme = fullfile(folder,'README.txt');
websave(file_readme,url_readme)

%% Scan Readme to obtain image types and details

% read Readme.txt
folder = fullfile(fileparts(which('xplor')),'demo','nasa_neo');
file_readme = fullfile(folder,'README.txt');
readme = fn_readtext(file_readme);

% process: separate lines, skip headers
% readme = strsplit(readme,'\n');
nheader = fn_find(@(str)strfind(str,'-------------'),readme,'first');
readme(1:nheader) = [];

% retrieve types and descriptions
ntype = length(readme);
alltypes = cell(1,ntype);
alltypes_descriptions = cell(1,ntype);
for i = 1:ntype
    [alltypes{i}, alltypes_descriptions{i}] = fn_regexptokens(readme{i},'^([^ ]*) +(.*[^ ]) *$');
end

%% Download full datasets! 
folder = fullfile(fileparts(which('xplor')),'demo','nasa_neo');
fn_mkdir(folder);

base_url = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/';

% start with smaller data set
types = { ...
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

% download loop
for type = types
    type_str = type{1}; % Get type as a string
    subfolder = fullfile(folder,type_str);
    fn_mkdir(subfolder) % Create folder in current Matlab path

    % list from web folder
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

%% Download specific files 
folder = fullfile(fileparts(which('xplor')),'demo','nasa_neo');
fn_mkdir(folder);

base_url = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/';

types = {'MOD_LSTAD_M' 'MOD_LSTAN_M' 'MOD_LSTD_M' 'MOD_LSTN_M'};
years_range = 2017:2019;

% download loop
for type = types
    type_str = type{1}; % Get type as a string
    subfolder = fullfile(folder,type_str);
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


%% Read data

% types = {'MOD_LSTAD_M' 'MOD_LSTAN_M' 'MOD_LSTD_M' 'MOD_LSTN_M'};
types = {'AURA_OZONE_M' 'CERES_INSOL_M' 'CERES_LWFLUX_M'}; %  'AVHRR_SST_M'
ntype = length(types);

desc = fn_map(@(type)alltypes_descriptions{strcmp(type,alltypes)},types);

% get year range by scanning which first and last file of each folder
years_range = zeros(ntype,2);
for ktype = 1:ntype
    type_str = types{ktype}; % Get type as a string
    subfolder = fullfile(folder,type_str);
    d = dir(fullfile(subfolder,'*.TIFF'));
    years_range(ktype,1) = str2double(fn_regexptokens(d(1).name,[type_str '_(\d{4})']));
    years_range(ktype,2) = str2double(fn_regexptokens(d(end).name,[type_str '_(\d{4})']));
end
disp([char(types) repmat(': ',ntype,1) num2str(years_range(:,1)) repmat(' -> ',ntype,1) num2str(years_range(:,2))])
years_range = min(years_range(:,1)):max(years_range(:,2));
nyear = length(years_range);

% spatial binning
nx = 3600/5;
ny = 1800/5;

% array
clear data
data = NaN([nx ny 12 nyear ntype]);

fn_progress('reading file',ntype*12*nyear)
kfile = 0;
for ktype = 1:ntype    
    type_str = types{ktype}; % Get type as a string
    subfolder = fullfile(folder,type_str);
    for kyear = 1:nyear
        year = years_range(kyear);
        for month = 1:12
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
                data(:,:,month,kyear,ktype) = x;
            else
                %fprintf('missing data for type %s, month %.2i/%i\n',type_str,month,year);
            end
        end
    end
end



