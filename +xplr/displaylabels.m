classdef displaylabels < xplr.graphnode
% displaylabels

    properties (Access='private')
        % parent xplr.viewdisplay object and other external objects
        D
        graph
        % internal
        height     % unit: normalized to axis size
        rotheight  % unit: normalized to axis size
        h
        movingdim
        movingclone
        listeners = struct('slice',[]);
        prevorgSetpos     % last organization seen by setPositions
    end
    properties (SetObservable=true)
        doImmediateDisplay = true;
    end
    properties (Dependent, Access='private')
        ha
    end
        
    methods
        function L = displaylabels(D)
            % contructor displaylabels
            
            % parent xplr.viewdisplay object
            L.D = D;
            
            % create labels and position them
            L.graph = D.graph;
            createLabels(L,'global')
            getHeights(L)
            if ~isempty(D.layoutID), setPositions(L), end
            % watch in axes size (no need to take care of deleting this
            % listener when L is deleted, because there is no reason that L
            % be deleted without D, D.ha, and the listener itself being
            % deleted as well)
            fn_pixelsizelistener(D.ha,@(u,e)updateLabels(L,'axsiz'))
            % note that changes in D lead to display updates through direct
            % method calls rather than notifications
        end
        function ha = get.ha(L)
            ha = L.D.ha;
        end
    end
    
    % Create and position labels
    methods (Access='private')
        function createLabels(L,flag,dims)
            % create labels
            if nargin<2, flag = 'global'; end
            if nargin<3 || strcmp(flag,'global'), dims = 1:L.D.nd; end
            [dims, dimIDs] = L.D.slice.dimensionNumberAndID(dims);
            curheaders = L.D.zslice.header;
            switch flag
                case 'global'
                    deleteValid(L.h)
                    L.h = gobjects(1,L.D.nd);
                otherwise
                    error('invalid flag ''%s''',flag)
            end
            for i = 1:length(dims)
                d = dims(i); dimID = dimIDs(i);
                str = curheaders(d).label;
                if ~isempty(curheaders(d).unit), str = [str ' (' curheaders(d).unit ')']; end
                L.h(d) = text('string',['  ' str '  '],'parent',L.ha, ...
                    'margin',1, ...
                    'backgroundcolor',[1 1 1]*.95,'units','normalized', ...
                    'UserData',dimID,'buttondownfcn',@(u,e)labelClick(L,u), ...
                    'UIContextMenu',uicontextmenu(L.D.V.hf,'callback',@(m,e)L.D.dimensionContextMenu(m,dimID)));
            end
        end
        function changeLabel(L,d)
            % change properties that need to when data undergoes a 'chgdim'
            % change
            
            % label
            head = L.D.zslice.header(d);
            str = head.label;
            if ~isempty(head.unit), str = [str ' (' head.unit ')']; end
            set(L.h(d),'string',['  ' str '  '])
            set(L.h(d),'UserData',head.dimID)
        end
        function getHeights(L)
            % get heights
            fontsize = get(L.D.hp,'defaultuicontrolfontsize'); % unit: points
            hinch = 1.5*fontsize/72;
            hpix = hinch*get(0,'ScreenPixelsPerInch');
            axsiz = fn_pixelsize(L.ha);
            L.height = hpix/axsiz(2);     % normalized to axis size
            L.rotheight = hpix/axsiz(1);  % normalized to axis size
            L.prevorgSetpos = []; % next call to setPositions should have doloose set to false
        end
        function setPositions(L)
            persistent prevorg
            
            % current layout
            org = L.D.layout;
            
            % locations of dimension: note that some dimensions might not
            % be present in the layout, in particular when it was removed
            % from the layout in labelMove, but the change is not
            % definitive and the slice has not been updated yet
            dim_locations = L.D.layoutID.dim_locations;
            
            % visible labels and active dimensions
            sz = L.D.slice.sz; % slice size
            okdim = (sz>1) & ~fn_isemptyc(dim_locations);
            isactive = false(1,length(sz));
            isactive([L.D.activedim.x L.D.activedim.y]) = true;
            
            % do 'loose' update? (i.e. do not adjust position for
            % non-relevant coordinates)
            % NO MORE NEEDED, as automatic positionning should be good
            % enough now
            if isempty(L.movingdim)
                doloose = false; %isequal(org,prevorg);
                prevorg = org;
            else
                doloose = false;
                prevorg = [];
            end
            
            % steps in the direction orthogonal to the positioning one
            axis_pos = fn_pixelpos(L.D.ha);
            available_space = axis_pos(1:2)./axis_pos(3:4) - L.height;
            xvstep = min(1.5*L.height, available_space(2)/(1.5+length(org.x)));
            yhstep = min(1.5*L.rotheight, available_space(1)/(1.5+length(org.y)));
            
            % set positions 
            for d=1:length(sz)
                if ~isgraphics(L.h(d),'text')
                    % it can happen that labels do not exist yet when
                    % updateLabels(L,'axsiz') is invoked
                    prevorg = []; % positions will not be set, so do not store a 'prevlayout'!
                    continue
                end
                if ~okdim(d)
                    % do not display the label of singleton dimension
                    set(L.h(d),'visible','off')
                    if any(d==L.movingdim)
                        set(L.movingclone,'visible','off')
                    end
                else
                    f = dim_locations{d};
                    % fixed properties
                    set(L.h(d),'visible','on', ...
                        'rotation',fn_switch(f,{'y' 'yx'},90,0), ...
                        'horizontalalignment',fn_switch(f,{'x' 'y'},'center',{'xy' 'mergeddata'},'left','yx','right'), ...
                        'verticalalignment',fn_switch(f,'yx','bottom','middle'), ...
                        'EdgeColor',fn_switch(isactive(d),'k','none'))
                    % set position
                    switch f
                        case 'x'
                            i = find(org.x==d,1);
                            xpos = .5 + L.graph.labelPosition(d);
                            newpos = [xpos -(1.5+i)*xvstep];
                        case 'y'
                            i = find(org.y==d,1);
                            ypos = .5 + L.graph.labelPosition(d);
                            newpos = [-(1.5+i)*yhstep ypos];
                        case 'mergeddata'
                            i = find(org.mergeddata==d,1);
                            newpos = [-4*L.rotheight (length(org.mergeddata)+.5-i)*L.height];
                        case 'xy'
                            newpos = [0 1+1.5*L.height];
                        case 'yx'
                            newpos = [1+2.2*L.rotheight 1];
                    end
                    if doloose
                        pos = get(L.h(d),'position');
                        switch f
                            case 'x'
                                newpos(2) = pos(2); set(L.h(d),'position',newpos)
                            case 'y'
                                newpos(1) = pos(1); set(L.h(d),'position',newpos)
                            otherwise
                                set(L.h(d),'position',newpos)
                        end
                    else
                        if any(d==L.movingdim)
                            set(L.movingclone,'position',newpos, ...
                                'rotation',get(L.h(d),'rotation'), ...
                                'visible',true, ...
                                'horizontalalignment',get(L.h(d),'horizontalalignment'), ...
                                'verticalalignment',get(L.h(d),'verticalalignment'))
                        else
                            set(L.h(d),'position',newpos)
                        end
                    end
                end
            end
        end
    end
    
    % Update labels
    methods
        function updateLabels(L,flag,dim)
            if nargin<2, flag = 'pos'; end
            switch flag
                case 'global'
                    createLabels(L,'global')
                case 'chgdim'
                    % some specific properties need to be updated
                    for d=dim, changeLabel(L,d), end
                case 'axsiz'
                    if isempty(L.D.layoutID), return, end % happens sometimes at init because figure size changes for no clear reason
                    getHeights(L)
                case 'pos'
                case 'active'
                    % only mark labels as active or not
                    nd = L.D.nd;
                    isactive = false(1,nd);
                    isactive([L.D.activedim.x L.D.activedim.y]) = true;
                    for d=1:nd
                        set(L.h(d),'EdgeColor',fn_switch(isactive(d),'k','none'))
                    end
                    return
                otherwise
                    error('invalid flag ''%s''',flag)
            end
            setPositions(L)
        end
    end
    
    % Label click
    methods
        function labelClick(L,obj)
            % dimension number
            dimID = get(obj,'UserData');
            % which action
            switch get(L.D.V.hf,'SelectionType')
                case 'normal'
                    % move the label; if it is not moved, change active dim
                    labelMove(L,dimID)
            end
        end
        function labelMove(L,dimID,do_swap)
            % function labelMove(L,dimID,do_swap)
            %---
            % Input:
            % - dimID   identifier of the dimension of the label to move
            % - do_swap if this dimension is in x and is moved to y, allow
            %           swapping with elements that are in y (this results
            %           for example in transposing images)
            
            % input
            if nargin<3, do_swap = true; end
            
            % prepare for changing organization
            [prevlayoutID, layoutID_d, layoutID] = deal(L.D.layoutIDall); % previous, previous without d, current
            dlayout = layoutID.dim_location(dimID);
            didx = find(layoutID.(dlayout)==dimID);
            layoutID_d.(dlayout) = setdiff(layoutID.(dlayout),dimID,'stable');
            
            % only one dimension gan go to either xy or yx, and if image
            % display only one dimension can go to mergeddata=colordim
            % -> remember the current dimension in these location for a
            % swap
            xydimID = [layoutID_d.xy layoutID_d.yx];
            if strcmp(L.D.displaymode,'image') 
                colordimID = layoutID_d.mergeddata;
            else
                colordimID = []; 
            end
            
            % data
            sz = L.D.slice.sz; okdim = (sz>1);
            dim = L.D.slice.dimensionNumber(dimID);
            dimIDsok = [L.D.slice.header(okdim).dimID];
            
            % set thresholds: 
            % thresholds will first be computed while ignoring
            % singleton dimensions, then NaNs will be inserted at the
            % locations of these singleton dimensions
            % (x)
            xlayout = layoutID.x; xthr = .5+L.graph.labelPosition(xlayout); % first get threshold with d and singleton dimensions still included
            xthr(xlayout==dimID) = []; xlayout(xlayout==dimID) = [];        % remove dimension d
            okx = ismember(xlayout,dimIDsok); xthr(~okx) = NaN;             % make it impossible to insert to the right of a singleton or non-present dimension
            if do_swap && strcmp(dlayout,'y')  % refine intervals in case of swapping
                % not only x insertions, but also x/y swaps are possible
                % for example:
                %      pos =    0   x           x               1
                %   -> thr =       | |       |     |
                nokx = sum(okx);
                pos = [0 xthr(okx) 1]; % add absolute min and max
                xthrok = zeros(2,nokx);
                for ithr=1:nokx
                    dmax = min(pos(ithr+1)-pos(ithr),pos(ithr+2)-pos(ithr+1))/4; % maximal distance to an x-label to swap with this label
                    xthrok(:,ithr) = pos(ithr+1)+[-dmax dmax];
                end
                xthr = NaN(2,length(xlayout)); xthr(:,okx) = xthrok;
            end
            xthr = [-Inf row(xthr)]; 
            % (y)
            ylayout = layoutID.y; ythr = .5-L.graph.labelPosition(ylayout); % first get threshold with d and singleton dimensions still included
            ythr(ylayout==dimID) = []; ylayout(ylayout==dimID) = [];        % remove dimension d
            oky = ismember(ylayout,dimIDsok); ythr(~oky) = NaN;             % make it impossible to insert to the right of a singleton or non-present dimension
            if do_swap && strcmp(dlayout,'x')
                % not only y insertions, but also x/y swaps are possible
                noky = sum(oky);
                pos = [0 ythr(oky) 1];
                ythrok = zeros(2,noky);
                for ithr=1:noky
                    dmax = min(pos(ithr+1)-pos(ithr),pos(ithr+2)-pos(ithr+1))/4; % maximal distance to an x-label to swap with this label
                    ythrok(:,ithr) = pos(ithr+1)+[-dmax dmax];
                end
                ythr = NaN(2,length(ylayout)); ythr(:,oky) = ythrok;
            end
            ythr = [-Inf row(ythr)];
                        
            % moving out of the graph to the "filter" area
            view_control = L.D.V.C; % what a nice piece of code isn't it!?
            do_filter = false;
            
            % label object, make sure it will not be covered by data display
            L.movingdim = L.D.slice.dimensionNumber(dimID);
            obj = L.h(L.movingdim);
            uistack(obj,'top')
            
            % move
            L.movingclone = copyobj(obj,L.ha);
            set(L.movingclone,'color',[1 1 1]*.6,'edgecolor','none','BackgroundColor','none')
            uistack(L.movingclone,'bottom')
            % (execute movelabel once: this is needed when L.labelMove is
            % called from xplr.viewcontrol.move_filtered_dimension and
            % mouse cursor is actually at a different position than the
            % label)
            movelabel
            % (execute movelabel upon mouse motion)
            moved = fn_buttonmotion(@movelabel,L.D.V.hf,'moved?','pointer','hand');

            function movelabel
                p = get(L.ha,'currentpoint'); p = p(1,1:2)';
                pfig = fn_coordinates(L.ha,'a2p',p,'position'); % convert to 'normalized' unit in parent panel
                p = fn_coordinates(L.ha,'a2c',p,'position');    % convert to 'normalized' unit in axes
                % move object
                set(obj,'position',p)
                % special: going outside of graph, in the controls area
                % -> filter dimension, new layout has no more dimension dimID
                if pfig(1) < 0
                    if ~do_filter
                        do_filter = true;
                        view_control.show_inoperant_filter(dimID)
                    end
                else
                    if do_filter
                        do_filter = false;
                        view_control.remove_inoperant_filter(dimID)
                    end
                end
                % immediate display update
                immediate_display = L.doImmediateDisplay;
                % update organization and object location
                newlayoutID = layoutID_d;
                if pfig(1) < -2
                    % nothing to do: newlayoutID is already the current
                    % layout without dimension d
                    % however we cannot perform a complete display update
                    % because the slice is not recomputed yet
                    immediate_display = false;
                elseif p(2)<=0 && p(2)<=p(1)
                    % insert in x
                    idx = find(p(1)>=xthr,1,'last'); % never empty thanks to the -Inf
                    % x/y swap rather than insert?
                    if do_swap && strcmp(dlayout,'y')
                        is_swap = ~mod(idx,2);
                        idx = ceil(idx/2);
                    else
                        is_swap = false;
                    end
                    if is_swap
                        newlayoutID.x = [layoutID_d.x(1:idx-1) dimID layoutID_d.x(idx+1:end)];
                        newlayoutID.y = [layoutID_d.y(1:didx-1) layoutID_d.x(idx) layoutID_d.y(didx:end)];
                    else
                        newlayoutID.x = [layoutID_d.x(1:idx-1) dimID layoutID_d.x(idx:end)];
                    end
                elseif p(1)<=0 && p(2)<=.25
                    % goes in mergeddata
                    newlayoutID.mergeddata(end+1) = dimID;
                    if colordimID
                        % move away dimension that was occupying mergeddata
                        % location
                        tmp = newlayoutID.(dlayout);
                        newlayoutID.(dlayout) = [tmp(1:didx-1) colordimID tmp(didx:end)];
                    end
                elseif p(1)<=0
                    % insert in y
                    idx = find(1-p(2)>=ythr,1,'last'); % never empty thanks to the +Inf
                    % x/y swap rather than insert?
                    if do_swap && strcmp(dlayout,'x')
                        is_swap = ~mod(idx,2);
                        idx = ceil(idx/2);
                    else
                        is_swap = false;
                    end
                    if is_swap
                        newlayoutID.x = [layoutID_d.x(1:didx-1) layoutID_d.y(idx) layoutID_d.x(didx:end)];
                        newlayoutID.y = [layoutID_d.y(1:idx-1) dimID layoutID_d.y(idx+1:end)];
                    else
                        newlayoutID.y = [newlayoutID.y(1:idx-1) dimID newlayoutID.y(idx:end)];
                    end
                elseif any(p>=1)
                    % xy and yx
                    if p(2)>=p(1)
                        % (zone above the graph)
                        [newlayoutID.xy, newlayoutID.yx] = deal(dimID,[]);
                    else
                        % (zone to the right of the graph)
                        [newlayoutID.xy, newlayoutID.yx] = deal([],dimID);
                    end
                    if xydimID
                        % move away dimension that was occupying xy or yx
                        % (swap with dimID previous location)
                        tmp = newlayoutID.(dlayout);
                        newlayoutID.(dlayout) = [tmp(1:didx-1) xydimID tmp(didx:end)];
                    end
                else
                    % no change in layout, return
                    % note that there is no more zone for yx!, meaning that
                    % this display option is no longer accessible, it
                    % seemed to be too useless
                    return
                end
                % update organization (-> will trigger automatic display
                % update)
                if ~isequal(newlayoutID,layoutID)
                    layoutID = newlayoutID;
                    L.D.setLayoutID(layoutID,immediate_display)
                end
                drawnow update
            end
            
            % finish 
            L.movingdim = [];
            delete(L.movingclone)
            L.movingclone = [];
            if ~moved
                % label was not moved, make d the active x or y dimension
                % if appropriate
                makeDimActive(L.D,dimID,'toggle')
            elseif do_filter
                % apply the filter! -> this will cause a reslice and  a
                % global change in dimensions
                view_control.activate_inoperant_filter(dimID)
            elseif isequal(layoutID,prevlayoutID)
                % update label positions once more to put the one for dimension
                % d in place, but keep the "non-so-meaningful" coordinate
                % at its new position (maybe user moved the label to a
                % better place where it does not hide other information)
                setPositions(L)
            elseif L.doImmediateDisplay || isequal(layoutID,prevlayoutID)
                % update label positions once more to put the one for dimension
                % d in place
                setPositions(L)
            else
                % change organization (which will trigger data display
                % update)only now
                L.D.setLayoutID(layoutID,true); % second argument (true) is needed to force display update even though D.layout is already set to curlayout
            end
        end
    end
    
end


