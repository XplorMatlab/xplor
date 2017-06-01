function test(varargin)
% Example:
% test [same as 'test all']
% test xdata
% test filter slicer

flags = varargin;
if nargin==0
    flags = {'all'};
end

% execute all tests
for i=1:length(flags)
    disp(['TEST ' upper(flags{i})])
    eval(['test_' flags{i}])
end

%---
function test_all

test_header
test_xdata
test_filterAndPoint
test_view

%---
function test_header %#ok<*DEFNU>

header = xplr.header({'x' 3},{'y' 2},{'time' 's' 5 0 1},{'cond' {'a' 'b' 'a' 'b'}});

%---
function test_xdata

header = xplr.header({'x' 3},{'y' 2},{'time' 's' 5 0 1},{'cond' {'a' 'b' 'a' 'b'}});
x = reshape(1:120,[3 2 5 4]);
data = xplr.xdata(x,header);
disp(data)

%---
function test_filterAndPoint

header = xplr.header({'x' 3},{'y' 2},{'time' 's' 5 0 1},{'cond' {'a' 'b' 'a' 'b'}});
x = reshape(1:120,[3 2 5 4]);
data = xplr.xdata(x,header);
F = xplr.filterAndPoint(header(1),'indices');
F.updateSelection('all',{1:2 3},'origin',{'titi' 'tata'})
slice = F.operation(data,1);


%---
function test_view

header = xplr.header({'x' 3},{'y' 2},{'time' 's' 5 0 1},{'cond' {'a' 'b' 'a' 'b'}});
x = reshape(1:120,[3 2 5 4]);
data = xplr.xdata(x,header);
Z = xplr.zoomfilter(header(3));
S = xplr.slicer(data,3,Z);
V = xplr.view(S.slice);

% %% OLD CODE
% 
% %% view
% 
% V = xplr.view(data);
% S = V.slicer;
% 
% %%
% F1 = F;
% S.addFilter(1,F1)
% 
% %%
% F2 = xplr.filterAndPoint(header(2),'indices');
% S.addFilter(2,F2)
% 
% %%
% F2.updateSelection('all',{[1 2] [3 4]})
% 
% %%
% F1.updateSelection('new',{1:3})
% 
% %% lists
% 
% fn_figure('LISTS')
% xplr.list(F1,'in',subplot(121));
% xplr.list(F2,'in',subplot(122));
% 
% % %% display
% % 
% % figure(1)
% % subplot(121)
% % H1 = hybriddisplay(data,[1 2]);
% % subplot(122)
% % H2 = hybriddisplay(data,3);
% % 
% % %%
% % 
% % clf
% % a = axes('pos',[.1 .1 .2 0],'tickdir','out','ticklength',[.02 .1]);
% % b = axes('pos',[.1 .1 .8 0],'tickdir','out');
