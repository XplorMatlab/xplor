
%% Load data
if ~exist('dataall','var')
    load '../demo data/unemployment'
    % reduce to origin = TOTAL -> make the data 4D
    clear origin
    dataall = squeeze(dataall(:,:,10,:,:));
    % keep only male/female
    idx = [1 2];
    sex = sex(idx);
    dataall = dataall(idx,:,:,:);
    % keep only ages by 5-years periods
    idx = [1 7 9 15 16 17 20 21 25 27 28 30];
    age = age(idx);
    dataall = dataall(:,idx,:,:);
    % keep only countries (discard groups of countries)
    idx = setdiff(1:39,[9:11 15:17]);
    country = country(idx);
    dataall = dataall(:,:,idx,:);
    % normalize the data between -.5 and .5
    M = max(dataall(:)); m = min(dataall(:));
    dataall = (dataall-m)/(M-m);
end

%% Make the xdata object

header = xplrlight.header({'date' 'year' length(date) date(1) mean(diff(date))}, ...
    {'sex' sex},{'age' age},{'country' country});
data = xplrlight.xdata(permute(dataall,[4 1 2 3]),header);

%% 'fake' view object
hf = fn_figure('XPLOR LIGHT','color','w');
V = struct('hf',hf);
V.panels.display = uipanel('parent',hf);
V.slicer.slice = data;

%% Display!
D = xplrlight.viewdisplay(V)

%% Some notes

% The following comments mark specific code removal
%L1     labels
%L2     clipping control
%L3     events and listeners
%L4     image display
%L5     changes other than 'global'
%L6     "dimension pairing"
%L7     line color
%L      all other features, including zoom, mouse navigation