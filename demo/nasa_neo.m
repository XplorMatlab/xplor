
%% Scan readme to obtain image types and details

% download Readme
url_readme = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/README.txt';
readme = webread(url_readme);

% process: separate lines, skip headers
readme = strsplit(readme,'\n');
nheader = fn_find(@(str)strfind(str,'-------------'),readme,'first');
readme(1:nheader) = [];

% retrieve types and descriptions
ntype = length(readme);
alltypes = cell(1,ntype);
alltypes_descriptions = cell(1,ntype);
for i = 1:ntype
    [alltypes{i}, alltypes_descriptions{i}] = fn_regexptokens(readme{i},'^([^ ]*) +(.*[^ ]) *$');
end

%% Download files 
folder = fullfile(fileparts(which('xplor')),'demo','nasa_neo');
fn_mkdir(folder);

base_url = 'https://neo.sci.gsfc.nasa.gov/archive/geotiff.float/';

% start with smaller data set
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

types = {'MOD_LSTAD_M' 'MOD_LSTAN_M' 'MOD_LSTD_M' 'MOD_LSTN_M'};
desc = fn_map(@(type)alltypes_descriptions{strcmp(type,alltypes)},types);
years_range = 2017:2019;

ntype = length(types);
nyear = length(years_range);

% spatial binning
xbin = 2;
nx = 3600/xbin;
ny = 1800/xbin;

% array
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
                x = fn_bin(x,xbin);
                % put in array
                data(:,:,month,kyear,ktype) = x;
            else
                fprintf('missing data for type %s, month %.2i/%i\n',type_str,month,year);
            end
        end
    end
end



