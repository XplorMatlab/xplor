function test(varargin)
% Example:
% test [same as 'test all']
% test xdata
% test filter slicer

flags = varargin;
if nargin == 0
    flags = {'all'};
end

% execute all tests
for i=1:length(flags)
    disp(['TEST ', upper(flags{i})])
    eval(['test_', flags{i}])
end

%---
function test_private
% test th
strict_size(2, 3)

%---
function test_all

test_header
test_xdata
test_filter_and_point
test_view

%---
function test_header %#ok<*DEFNU>

header = xplr.Header({'x', 3}, {'y', 2}, {'time', 's', 5, 1, 0}, {'cond', {'a',
 'b', 'a', 'b'}});

%---
function test_xdata

header = xplr.Header({'x', 3}, {'y', 2}, {'time', 's', 5, 1, 0}, {'cond', {'a', 'b', 'a', 'b'}});
x = reshape(1:120, [3, 2, 5, 4]);
data = xplr.XData(x, header);
disp(data)

%---
function test_filter_and_point

header = xplr.Header({'x', 3}, {'y', 2}, {'time', 's', 5, 1, 0}, {'cond', {'a', 'b', 'a', 'b'}});
x = reshape(1:120, [3, 2, 5, 4]);
data = xplr.XData(x, header);
F = xplr.FilterAndPoint(header(1));
F.update_selection('all', {1:2, 3}, 'origin', {'titi', 'tata'})
slice = F.operation(data, 1);

%---
function test_view

header = xplr.Header({'x', 3}, {'y'  2}, {'time', 's', 5, 1, 0}, {'cond', {'a', 'b', 'a', 'b'}});
x = reshape(1:120, [3, 2, 5, 4]);
data = xplr.XData(x, header);
V = xplr.View(data);

%---
function test_zoom_central
%%

% two displays with y coordinate not sampled the same way
load 'clown';
x1 = xplr.XData(X, {{'x', 'px'}, {'y', 'cm', 2/320, -1+1/320}});
x2 = xplr.XData(X(:, 1:2:end), {{'x', 'px'}, {'y', 'cm', 2/160, -1+1/160}});

V1 = xplor(x1);
V2 = xplor(x2);

error 'much of the code below is not working any more due to changes in code'
% % connect the 2 zoom filters in the y coordinate
% ZF1 = V1.D.zoomslicer.getFilterAndPoint(2);
% ZF2 = V2.D.zoomslicer.getFilterAndPoint(2);
% 
% %%
% ZC = xplr.zoomcentral('y', 'px');
% 
% % notify(ZF1,'changed_zoom',xplr.eventinfo('zoom',0))
% 
% ZC.connectZoomFilter(ZF1);
% ZC.connectZoomFilter(ZF2);


% %% OLD CODE
% 
% %% view
% 
% V = xplr.view(data);
% S = V.slicer;
% 
% %%
% F1 = F;
% S.add_filter(1,F1)
% 
% %%
% F2 = xplr.filterAndPoint(header(2));
% S.add_filter(2,F2)
% 
% %%
% F2.update_selection('all',{[1 2] [3 4]})
% 
% %%
% F1.update_selection('new',{1:3})
% 
% %% lists
% 
% brick.figure('LISTS')
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
% % a = axes('position',[.1 .1 .2 0],'tickdir','out','ticklength',[.02 .1]);
% % b = axes('position',[.1 .1 .8 0],'tickdir','out');
