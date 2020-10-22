function d = demo_data_folder

if isdeployed
    d  = brick.userconfig('xplor_demo_data_folder');
    if ~exist(d, 'dir')
        d = brick.getdir('Please select the XPLOR demo data folder');
        brick.userconfig('xplor_demo_data_folder', d)
    end
else
    d = fullfile(fileparts(which('xplor')), 'demo data');
end
