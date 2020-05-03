function test_selections
%%
V = evalin('base','V');

D = V.D;
N = D.navigation;
N.selectiondim = [1 2];

x = V.slice.data(:,:,1);

X = fourd(x,'2d');

add_listener(X.SI,'ChangeView',@changeview);

    function changeview(u,e)

    % we are interested only in 'selection' changes
    if ~strcmp(e.flag, 'selection'), return, end

    % set navigation selection and display it
    t = X.SI.selection.t;
    if isempty(t)
        value = [];
    else
        value = X.SI.selection.t.set;
    end
    
    N.selection = value;
    N.displayselection()
    
    %     % for later...
    %     N.updateselection(e.selflag, e.ind, e.value)
    
    end

end