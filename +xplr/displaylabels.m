classdef displaylabels < xplr.graphnode
% displaylabels

    properties (Access='private')
        % parent xplr.viewdisplay object and other external objects
        D
        graph
        % internal
        height
        rotheight
        h
        movingdim
        listeners = struct('slice',[]);
        prevorgSetpos     % last organization seen by setPositions
    end
    properties
        doImmediateDisplay = false;
    end
    properties (Dependent, Access='private')
        ha
        org
    end
        
    methods
        function L = displaylabels(D)
            % contructor displaylabels
            
            % parent xplr.viewdisplay object
            L.D = D;
            
            L.graph = D.graph;
            createLabels(L,'global')
            getHeights(L)
            if ~isempty(D.org), setPositions(L), end
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
        function org = get.org(L)
            org = L.D.org;
        end
    end
    
    % Create and position labels
    methods (Access='private')
        function createLabels(L,flag,dim)
            % create labels
            if nargin<2, flag = 'global'; end
            if nargin<3 || strcmp(flag,'global'), dim = 1:L.D.nd; end
            curheaders = L.D.zslice.header;
            switch flag
                case 'global'
                    deleteValid(L.h)
                    L.h = gobjects(1,L.D.nd);
                case 'insertdim'
                    if ~all(diff(dim)==1), error 'not implemented yet', end
                    L.h = [L.h(1:dim(1)-1) gobjects(1,length(dim)) L.h(dim(1):end)];
                case 'rmdim'
                    deleteValid(L.h(dim))
                    L.h(dim) = [];
                    return
                case 'permdim'
                    L.h = L.h(dim);
                    return
                otherwise
                    error('invalid flag ''%s''',flag)
            end
            for d=dim
                str = curheaders(d).label;
                if curheaders(d).ismeasure, str = [str ' (' curheaders(d).unit ')']; end
                L.h(d) = text('string',['  ' str '  '],'parent',L.ha, ...
                    'backgroundcolor',[1 1 1]*.95,'units','normalized', ...
                    'userdata',d,'buttondownfcn',@(u,e)labelClick(L,u));
            end
        end
        function changeLabel(L,d)
            % change properties that need to when data undergoes a 'chgdim'
            % change
            
            % label
            head = L.D.zslice.header(d);
            str = head.label;
            if head.ismeasure, str = [str ' (' head.unit ')']; end
            set(L.h(d),'string',['  ' str '  '])
        end
        function getHeights(L)
            % get heights
            fontsize = get(L.D.hp,'defaultuicontrolfontsize'); % unit: points
            hinch = 1.5*fontsize/72;
            hpix = hinch*get(0,'ScreenPixelsPerInch');
            axsiz = fn_pixelsize(L.ha);
            L.height = hpix/axsiz(2);
            L.rotheight = hpix/axsiz(1);
            L.prevorgSetpos = []; % next call to setPositions should have doloose set to false
        end
        function setPositions(L)
            % do 'loose' update? (i.e. do not adjust position for
            % non-relevant coordinates)
            curorg = L.org;
            if isempty(L.movingdim)
                doloose = isequal(curorg,L.prevorgSetpos);
                L.prevorgSetpos = curorg;
            else
                doloose = false;
                L.prevorgSetpos = [];
            end
            
            % do not display labels of singleton dimensions
            sz = L.D.slice.sz; % slice size
            okdim = (sz>1);
            %if ~isempty(curorg.x), okdim(curorg.x(1)) = true; end % first x dimension must be visible regardless of its being singleton or not
            %if strcmp(L.D.displaymode,'image') && ~isempty(curorg.y), okdim(curorg.y(1)) = true; end % similar comment
            xorgok = curorg.x(okdim(curorg.x));
            yorgok = curorg.y(okdim(curorg.y));
            ystatok = curorg.ystatic(okdim(curorg.ystatic));
            orgperdim = org2perdim(curorg,length(sz));
            isactive = false(1,length(sz)); 
            isactive([L.D.activedim.x L.D.activedim.y]) = true;
            
            % set positions 
            for d=1:length(sz)
                if ~isgraphics(L.h(d),'text')
                    % it can happen that labels do not exist yet when
                    % updateLabels(L,'axsiz') is invoked
                    L.prevorgSetpos = []; % positions will not be set, so do not store a 'prevorg'!
                    continue
                end
                if ~okdim(d)
                    % do not display the label of singleton dimension
                    set(L.h(d),'visible','off')
                else
                    f = orgperdim{d};
                    % fixed properties
                    set(L.h(d),'visible','on', ...
                        'rotation',fn_switch(f,{'y' 'yx'},90,0), ...
                        'horizontalalignment',fn_switch(f,{'x' 'y'},'center',{'xy' 'ystatic'},'left','yx','right'), ...
                        'verticalalignment',fn_switch(f,{'y' 'yx'},'bottom','middle'), ...
                        'EdgeColor',fn_switch(isactive(d),'k','none'))
                    % set position
                    switch f
                        case 'x'
                            i = find(xorgok==d,1);
                            xpos = .5 + L.graph.labelPosition(d);
                            newpos = [xpos -(1+i)*L.height];
                        case 'y'
                            i = find(yorgok==d,1);
                            ypos = .5 + L.graph.labelPosition(d);
                            newpos = [-(1+i)*L.rotheight ypos];
                        case 'ystatic'
                            i = find(ystatok==d,1);
                            newpos = [-4*L.rotheight (length(curorg.ystatic)+.5-i)*L.height];
                        case 'xy'
                            newpos = [0 1+1.5*L.height];
                        case 'yx'
                            newpos = [-(length(yorgok)+2)*L.rotheight 1];
                    end
                    if doloose
                        pos = get(L.h(d),'pos');
                        switch f
                            case 'x'
                                newpos(2) = pos(2); set(L.h(d),'pos',newpos)
                            case 'y'
                                newpos(1) = pos(1); set(L.h(d),'pos',newpos)
                            otherwise
                                set(L.h(d),'pos',newpos)
                        end
                    elseif ~any(d==L.movingdim)
                        set(L.h(d),'pos',newpos)
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
                case {'insertdim' 'rmdim' 'permdim'}
                    createLabels(L,flag,dim)
                case 'chgdim'
                    % some specific properties need to be updated
                    for d=dim, changeLabel(L,d), end
                case 'axsiz'
                    if isempty(L.org), return, end % happens sometimes at init because figure size changes for no clear reason
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
            d = get(obj,'UserData');
            % which action
            switch get(L.D.V.hf,'SelectionType')
                case 'normal'
                    % move the label; if it is not moved, change active dim
                    labelMove(L,d,obj)
                case 'alt'
                    % show context menu
                    L.D.V.context.raise('label',d)
            end
        end
        function labelMove(L,d,obj)
            % prepare for changing organization
            [prevorg curorg org0] = deal(L.org);
            F = fieldnames(L.org);
            for i=1:length(F)
                f = F{i};
                if ismember(d,L.org.(f))
                    dorg = f; didx = find(L.org.(f)==d);
                    org0.(f) = setdiff(L.org.(f),d,'stable');
                    break
                elseif i==length(F)
                    warning('dimension %i not found in current organization',d)
                    dorg = '';
                end
            end
            % (if d is assigned to xy or yx and there is already a
            % dimension there, this one will need to be moved away)
            xydim = [org0.xy org0.yx]; 
            % prepare the global variables used in 'setthresholds' and
            % 'movelabel'
            sz = L.D.slice.sz; okdim = (sz>1);
            xthr = []; ythr = [];
            setthresholds()
            function setthresholds
                % thresholds will first be computed while ignoring
                % singleton dimensions, then NaNs will be inserted at the
                % locations of these singleton dimensions
                % (x)
                xorg = curorg.x; xthr = .5+L.graph.labelPosition(xorg); % first get threshold with d and singleton dimensions still included
                xthr(xorg==d) = []; xorg(xorg==d) = [];                 % remove dimension d
                okx = okdim(xorg); xthr(~okx) = NaN;                    % make it impossible to insert to the right of a singleton dimension
                if strcmp(dorg,'y')                                     % refine intervals in case of swapping
                    % not only x insertions, but also x/y swaps are possible
                    nokx = sum(okx);
                    tmp = [0 xthr(okx) 1];
                    xthrok = zeros(2,nokx);
                    for ithr=1:nokx
                        dmax = min(tmp(ithr+1)-tmp(ithr),tmp(ithr+2)-tmp(ithr+1))/4; % maximal distance to an x-label to swap with this label
                        xthrok(:,ithr) = tmp(ithr+1)+[-dmax dmax];
                    end
                    xthr = NaN(2,length(xorg)); xthr(:,okx) = xthrok;
                end
                xthr = [-Inf row(xthr)];
                % (y)
                yorg = curorg.y; ythr = .5+L.graph.labelPosition(yorg); % first get threshold with d and singleton dimensions still included
                ythr(yorg==d) = []; yorg(yorg==d) = [];                 % remove dimension d
                oky = okdim(yorg); ythr(~oky) = NaN;                    % make it impossible to insert to the right of a singleton dimension
                if strcmp(dorg,'x')
                    % not only y insertions, but also x/y swaps are possible
                    noky = sum(oky);
                    tmp = [0 ythr(oky) 1];
                    ythrok = zeros(2,noky);
                    for ithr=1:noky
                        dmax = min(tmp(ithr)-tmp(ithr+1),tmp(ithr+1)-tmp(ithr+2))/4; % maximal distance to an x-label to swap with this label
                        ythrok(:,ithr) = tmp(ithr+1)+[dmax -dmax];
                    end
                end
                ythr = [+Inf row(ythr)];
            end
            % make sure label will not be covered by data display
            uistack(obj,'top')
            % move
            L.movingdim = d;
            %movelabel()
            moved = fn_buttonmotion(@movelabel,L.D.V.hf,'moved?');
            function movelabel
                p = get(L.ha,'currentpoint');
                p = fn_coordinates(L.ha,'a2c',p(1,1:2)','position'); % convert to 'normalized' unit
                % move object
                set(obj,'pos',p)
                % update organization and object location
                neworg = org0;
                if p(2)<=0 && p(2)<=p(1)
                    % insert in x
                    idx = find(p(1)>=xthr,1,'last'); % never empty thanks to the -Inf
                    % x/y swap rather than insert?
                    if strcmp(dorg,'y')
                        doswap = ~mod(idx,2);
                        idx = ceil(idx/2);
                    else
                        doswap = false;
                    end
                    if doswap
                        neworg.x = [org0.x(1:idx-1) d org0.x(idx+1:end)];
                        neworg.y = [org0.y(1:didx-1) org0.x(idx) org0.y(didx:end)];
                    else
                        neworg.x = [org0.x(1:idx-1) d org0.x(idx:end)];
                    end
                    %while ~okdim(neworg.x(1)), neworg.x = neworg.x([2:end 1]); end % do not let the first dimension become singleton (note that d is non-singleton, so the loop will not be infinite)
                elseif p(1)<=0 && strcmp(L.D.displaymode,'time courses') && p(2)<=.25
                    % goes in ystatic
                    neworg.ystatic(end+1) = d;
                elseif p(1)<=0
                    % insert in y
                    idx = find(p(2)<=ythr,1,'last'); % never empty thanks to the +Inf
                    % x/y swap rather than insert?
                    if strcmp(dorg,'x')
                        doswap = ~mod(idx,2);
                        idx = ceil(idx/2);
                    else
                        doswap = false;
                    end
                    if doswap
                        neworg.x = [org0.x(1:didx-1) org0.y(idx) org0.x(didx:end)];
                        neworg.y = [org0.y(1:idx-1) d org0.y(idx+1:end)];
                    else
                        neworg.y = [neworg.y(1:idx-1) d neworg.y(idx:end)];
                    end
                    %while ~okdim(neworg.y(1)), neworg.y = neworg.y([2:end 1]); end % do not let the first dimension become singleton (note that d is non-singleton, so the loop will not be infinite)
                else
                    % xy or yx
                    if p(1)>=p(2) || p(1)>=(1-p(2))
                        % xy arrangement is more usual than yx, therefore
                        % give it more 'space
                        neworg.xy = d; neworg.yx = [];
                    else
                        neworg.xy = []; neworg.yx = d;
                    end
                    if xydim
                        % move away dimension that was occupying xy or yx
                        tmp = neworg.(dorg);
                        neworg.(dorg) = [tmp(1:didx-1) xydim tmp(didx:end)];
                    end
                end
                % update organization (-> will trigger automatic display
                % update)
                if ~isequal(neworg,curorg)
                    curorg = neworg;
                    L.D.setOrg(curorg,L.doImmediateDisplay)
                end
                drawnow update
                % (do not) update position thresholds (this was nice before
                % the swap was introduced, now it make the whole think too
                % shaky)
                % setthresholds
            end
            % finish 
            L.movingdim = [];
            if ~moved
                % label was not moved, make d the active x or y dimension
                % if appropriate
                makeDimActive(L.D,d,'toggle')
            elseif isequal(curorg,prevorg)
                % update label positions once more to put the one for dimension
                % d in place, but keep the "non-so-meaningful" coordinate
                % at its new position (maybe user moved the label to a
                % better place where it does not hide other information)
                setPositions(L)
            elseif L.doImmediateDisplay || isequal(curorg,prevorg)
                % update label positions once more to put the one for dimension
                % d in place
                setPositions(L)
            else
                % change organization (which will trigger data display
                % update)only now
                L.D.setOrg(curorg,true); % second argument (true) is needed to force display update even though D.org is already set to curorg
            end
        end
    end
    
end


%---
function orgperdim = org2perdim(org,nd)

orgperdim = cell(1,nd);
F = fieldnames(org);
for i=1:length(F)
    f = F{i};
    orgperdim(org.(f)) = {f};
end

end


