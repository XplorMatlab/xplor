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
    properties
        doImmediateDisplay = false;
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
            dim_locations = L.D.layoutID.dim_locations;
            
            % visible labels and active dimensions
            sz = L.D.slice.sz; % slice size
            okdim = (sz>1);
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
                else
                    f = dim_locations{d};
                    % fixed properties
                    set(L.h(d),'visible','on', ...
                        'rotation',fn_switch(f,{'y' 'yx'},90,0), ...
                        'horizontalalignment',fn_switch(f,{'x' 'y'},'center',{'xy' 'ystatic'},'left','yx','right'), ...
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
                        case 'ystatic'
                            i = find(org.ystatic==d,1);
                            newpos = [-4*L.rotheight (length(org.ystatic)+.5-i)*L.height];
                        case 'xy'
                            newpos = [0 1+1.5*L.height];
                        case 'yx'
                            newpos = [-(length(org.y)+2)*L.rotheight 1];
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
                    labelMove(L,dimID,obj)
            end
        end
        function labelMove(L,dimID,obj)
            % prepare for changing organization
            [prevlayoutID, layoutID_d, layoutID] = deal(L.D.layoutIDmemory); % previous, previous without d, current
            dlayout = layoutID.dim_location(dimID);
            didx = find(layoutID.(dlayout)==dimID);
            layoutID_d.(dlayout) = setdiff(layoutID.(dlayout),dimID,'stable');
            % (if d is assigned to xy or yx and there is already a
            % dimension there, this one will need to be moved away)
            xydim = [layoutID_d.xy layoutID_d.yx];
            
            % data
            sz = L.D.slice.sz; okdim = (sz>1);
            dimIDsok = [L.D.slice.header(okdim).dimID];
            
            % set thresholds: 
            % thresholds will first be computed while ignoring
            % singleton dimensions, then NaNs will be inserted at the
            % locations of these singleton dimensions
            % (x)
            xlayout = layoutID.x; xthr = .5+L.graph.labelPosition(xlayout); % first get threshold with d and singleton dimensions still included
            xthr(xlayout==dimID) = []; xlayout(xlayout==dimID) = [];        % remove dimension d
            okx = ismember(xlayout,dimIDsok); xthr(~okx) = NaN;        % make it impossible to insert to the right of a singleton or non-present dimension
            if strcmp(dlayout,'y')                                     % refine intervals in case of swapping
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
            if strcmp(dlayout,'x')
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
            
            % make sure label will not be covered by data display
            uistack(obj,'top')
            % move
            L.movingdim = dimID;
            L.movingclone = copyobj(obj,L.ha);
            set(L.movingclone,'color',[1 1 1]*.6,'edgecolor','none','BackgroundColor','none')
            uistack(L.movingclone,'bottom')
            %movelabel()
            moved = fn_buttonmotion(@movelabel,L.D.V.hf,'moved?');
            
            function movelabel
                p = get(L.ha,'currentpoint');
                p = fn_coordinates(L.ha,'a2c',p(1,1:2)','position'); % convert to 'normalized' unit
                % move object
                set(obj,'position',p)
                % update organization and object location
                newlayoutID = layoutID_d;
                if p(2)<=0 && p(2)<=p(1)
                    % insert in x
                    idx = find(p(1)>=xthr,1,'last'); % never empty thanks to the -Inf
                    % x/y swap rather than insert?
                    if strcmp(dlayout,'y')
                        doswap = ~mod(idx,2);
                        idx = ceil(idx/2);
                    else
                        doswap = false;
                    end
                    if doswap
                        newlayoutID.x = [layoutID_d.x(1:idx-1) dimID layoutID_d.x(idx+1:end)];
                        newlayoutID.y = [layoutID_d.y(1:didx-1) layoutID_d.x(idx) layoutID_d.y(didx:end)];
                    else
                        newlayoutID.x = [layoutID_d.x(1:idx-1) dimID layoutID_d.x(idx:end)];
                    end
                elseif p(1)<=0 && strcmp(L.D.displaymode,'time courses') && p(2)<=.25
                    % goes in ystatic
                    newlayoutID.ystatic(end+1) = dimID;
                elseif p(1)<=0
                    % insert in y
                    idx = find(1-p(2)>=ythr,1,'last'); % never empty thanks to the +Inf
                    % x/y swap rather than insert?
                    if strcmp(dlayout,'x')
                        doswap = ~mod(idx,2);
                        idx = ceil(idx/2);
                    else
                        doswap = false;
                    end
                    if doswap
                        newlayoutID.x = [layoutID_d.x(1:didx-1) layoutID_d.y(idx) layoutID_d.x(didx:end)];
                        newlayoutID.y = [layoutID_d.y(1:idx-1) dimID layoutID_d.y(idx+1:end)];
                    else
                        newlayoutID.y = [newlayoutID.y(1:idx-1) dimID newlayoutID.y(idx:end)];
                    end
                else
                    % xy or yx
                    if p(1)>=p(2) || p(1)>=(1-p(2))
                        % xy arrangement is more usual than yx, therefore
                        % give it more 'space
                        newlayoutID.xy = dimID; newlayoutID.yx = [];
                    else
                        newlayoutID.xy = []; newlayoutID.yx = dimID;
                    end
                    if xydim
                        % move away dimension that was occupying xy or yx
                        tmp = newlayoutID.(dlayout);
                        newlayoutID.(dlayout) = [tmp(1:didx-1) xydim tmp(didx:end)];
                    end
                end
                % update organization (-> will trigger automatic display
                % update)
                if ~isequal(newlayoutID,layoutID)
                    layoutID = newlayoutID;
                    L.D.setLayoutMemory(layoutID,L.doImmediateDisplay)
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
                L.D.setLayoutMemory(layoutID,true); % second argument (true) is needed to force display update even though D.layout is already set to curlayout
            end
        end
    end
    
end


