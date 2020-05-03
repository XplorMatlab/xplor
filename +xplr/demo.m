function demo
% Type xplr.demo to access a choice of demos of the XPLOR toolbox


demo_list = {'Intrinsic imaging'};
demo_fun = {@demo_intrinsic};
% demo_idx = listdlg('Name','XPLOR demo','PromptString','Select demo', ...
%     'ListString',demo_list,'SelectionMode','single');
demo_idx = 1;
if isempty(demo_idx), return, end

demo_dir = fullfile(fileparts(which('xplor')), 'demo data');
cd(demo_dir)
feval(demo_fun{demo_idx})

%---
function demo_intrinsic

%%
load intrinsic

fn_figure('Intrinsic imaging - Main', [1023, 598], 'color', 'w')

axes('position',[.02, .52, .3, .46])
a = fourd(img, '2d', 'labels', {'x', 'y'}, 'units', {'um', 'um'}, 'mat', [20/6, 20/6]);
a.D.do_labels = false;

axes('position', [.02, .02, .3, .46])
b = fourd(x, '2d', ...
    'labels', {'x', 'y', 'time', 'condition', 'repetition'}, ...
    'units', {'um', 'um', 's', {'stimulation', 'control'}, ''}, ...
    'mat', {[20, 20, .2, 1, 1], [-20, -20, -.2, 0, 0]});
b.D.do_labels = false;
G = b.G;

axes('position', [.38, .1, .6, .86])
c = fourd(x, '2dplot', G, 'dimsplus', []);

X = explor(x, G);
