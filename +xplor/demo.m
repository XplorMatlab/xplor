function demo
% Type xplr.demo to access a choice of demos of the XPLOR toolbox

% List of available demos
available_demos = list_available_demos;
demo_list = available_demos(:, 1);
demo_fun = available_demos(:, 2);

% Execute one demo or display the full list
hf = brick.figure('XPLOR demo', [240 400], 'Menubar', 'none', 'Toolbar', 'none');
h_label = .05;
uicontrol('style', 'text', 'parent', hf, ...
    'units', 'normalized', 'position', [0 1-h_label 1 h_label], ...
    'string', 'Select demo to run', 'fontsize', 12, 'horizontalalignment', 'left');
uicontrol('style', 'listbox', 'parent', hf, ...
    'units', 'normalized', 'position', [0 .5 1 .5-h_label], ...
    'string', demo_list, ...
    'callback', @(u, e)do_demo(get(u, 'value')));
uicontrol('style', 'text', 'parent', hf, ...
    'units', 'normalized', 'position', [0 .5-h_label 1 h_label], ...
    'string', 'Show demo code', 'fontsize', 12, 'horizontalalignment', 'left');
uicontrol('style', 'listbox', 'parent', hf, ...
    'units', 'normalized', 'position', [0 0 1 .5-h_label], ...
    'string', demo_list, ...
    'callback', @(u, e)edit(char(demo_fun{get(u, 'value')})));

%---
function available_demos = list_available_demos

available_demos = {
    'XPLOR logo',           @xplor.demo.logo
    'Flow',                 @xplor.demo.flow
    'Astronomy',            @xplor.demo.hubble_pillars_of_the_creation
    'Coronavirus evolution',    @xplor.demo.coronavirus
    'Intrinsic imaging',    @xplor.demo.intrinsic
    '2-photon microscopy',  @xplor.demo.neuron_movies
    'Earth temperature',    @xplor.demo.earth_temperature
    };

%---
function do_demo(demo_idx)

if isempty(demo_idx), return, end

% Close all XPLOR windows
close(findall(0, 'type', 'figure', 'tag', 'XPLOR'))

% Do demo
available_demos = list_available_demos;
demo_fun = available_demos(:, 2);
demo_fun{demo_idx}()

