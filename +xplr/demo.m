function demo
% Type xplr.demo to access a choice of demos of the XPLOR toolbox


demolist = {'Intrinsic imaging'};
demofun = {@demo_intrinsic};
% demoidx = listdlg('Name','XPLOR demo','PromptString','Select demo', ...
%     'ListString',demolist,'SelectionMode','single');
demoidx = 1;
if isempty(demoidx), return, end

demodir = fullfile(fileparts(which('xplor')),'demo data');
cd(demodir)
feval(demofun{demoidx})

%---
function demo_intrinsic

%%
load intrinsic

fn_figure('Intrinsic imaging - Main',[1023 598],'color','w')

axes('pos',[.02 .52 .3 .46])
a = fourd(img,'2d','labels',{'x' 'y'},'units',{'um' 'um'},'mat',[20/6 20/6]);
a.D.dolabels = false;

axes('pos',[.02 .02 .3 .46])
b = fourd(x,'2d', ...
    'labels',{'x' 'y' 'time' 'condition' 'repetition'}, ...
    'units',{'um' 'um' 's' {'stimulation' 'control'} ''}, ...
    'mat',{[20 20 .2 1 1] [-20 -20 -.2 0 0]});
b.D.dolabels = false;
G = b.G;

axes('pos',[.38 .1 .6 .86])
c = fourd(x,'2dplot',G,'dimsplus',[]);

X = explor(x,G);

