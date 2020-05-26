close all hidden
clear all

%% XPLOR FOLDER

xplor_folder = fileparts(which('xplor'));

%% NEURONS

if eval('false')
    %% ()
    s = fn_loadvar(fullfile(xplor_folder,'demo data','neurons_movie.mat'));
    x = xplr.xdata(s.data, ...
        {{'x' 'um' s.pixel_size}, {'y' 'um' s.pixel_size}, {'time' 's' 1/s.frame_rate}, 'trial'}, ...
        'NEURONS');
    
    V1 = xplor(x);
    V1.hf.Position = [24   389   560   473];
    V1.C.dimaction('viewdim',{'x' 'y'})
    V1.D.colormap.cmapdef = 'gray';
    V1.D.navigation.crosscolor = 'white';
    V1.D.navigation.selectiondim = {'x' 'y'};
    
    V2 = xplor(x);
    V2.hf.Position = [591   370   685   645];
    V2.C.dimaction('viewdim','time')
    V2.D.navigation.selectiondim = 'time';
end

%% NEO

if eval('false')
    %% ()
    s = fn_loadvar(fullfile(xplor_folder,'demo data','NEO Earth Temperature.mat'));
    
    %%
    data = s.temperatures(1+s.temperatures_indices);
    x = xplr.xdata(data, ...
        {{'latitude' '°' diff(s.latitude(1:2)) s.latitude(1)}, ...
        {'longitude' '°' diff(s.longitude(1:2)) s.longitude(1)}, ...
        {'months' s.months}}, ...
        'TEMPERATURE');
    %%
    V1 = xplor(x);
    V1.hf.Position = [686   361   592   640];
    V1.C.dimaction('viewdim',{'latitude' 'longitude'})
    V1.D.navigation.selectiondim = {'latitude' 'longitude'};
    %%
    V2 = xplor(x);
    V2.hf.Position = [566   361   712   640];
    V2.C.dimaction('viewdim','months')
    V2.D.navigation.selectiondim = 'months';
    
end


