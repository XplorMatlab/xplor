function wizard

fig_siz = [250 80];
screen_pos = get(0, 'ScreenSize');
W.hf = uifigure('name', 'XPLOR', 'pos', [5 screen_pos(4)-fig_siz(2)-30 fig_siz], ...
    'menubar', 'none', 'resize', 'off');
h = uihtml(W.hf);
h.HTMLSource = fullfile(fileparts(which('xplor')), '+xplr/wizard.html');
h.DataChangedFcn = @(src,event)eval(h.Data);
h.Position = [1 1 fig_siz];

%---
function button_press(name)

if ~isempty(name)
    eval(name)
end

%---
function LoadCSVData(file)

if nargin < 1
    file = brick.getfile('*.csv', 'Select data files to XPLOR');
    if isequal(file, 0), return, end
end
data = io.read_table(file);
xplor(data)

%---
function ShowDemos()

xplor.demo
