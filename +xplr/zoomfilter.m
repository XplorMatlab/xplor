classdef zoomfilter < xplr.dataOperand
    % function Z = zoomfilter(headerin[,zoom[,bin]])
    %---
    % defines zooming, but also binning
   
    properties (SetAccess='private')
        % filter definition
        zoom = ':'      % ':' or [istart istop] with 1 <= istart <= istop <= headerin.n
        bin = 1
        % output
        indicesin       % data points that will be extracted
        indicesout      % data positions after zooming AND BINING (is equal to indicesin only if bin=1)
    end
    properties (Dependent, SetAccess='private')
        zoomvalue       % [istart istop]
    end
    
    % Setting and updating filter
    methods
        function Z = zoomfilter(headerin,zoom,bin)
            % invalid zoomfilter
            if nargout==0, return, end

            % input
            if nargin<3, zoom = ':'; end
            if nargin<4, bin = 1; end
            if ~isscalar(headerin)
                Z = xplr.zoomfilter.empty(1,0);
                for i=1:length(headerin)
                    Z(i) = xplr.zoomfilter(headerin(i));
                end
                return
            end
            
            % input header
            Z.headerin = headerin;
            
            % operation definition
            Z.zoom = zoom;
            Z.bin = bin;
            
            % set indices and output header
            prepareFilter(Z,true,true)
        end
        function setZoom(Z,zoom,bin)
            chgzoom = ~isequal(zoom,Z.zoom);
            chgbin = (nargin>=3 && bin~=Z.bin);
            if ~chgzoom && ~chgbin, return, end
            % check and assign zoom
            if strcmp(zoom,':')
                Z.zoom = zoom;
            else
                if ~isnumeric(zoom) || length(zoom)~=2, error 'wrong zoom value', end
                zoom = [max(.5,zoom(1)) min(Z.headerin.n+.5,zoom(2))];
                if diff(zoom)<=0, zoom = ':'; end % invalid zoom -> zoom reset
                Z.zoom = zoom;
            end
            % assign bin and update output
            if nargin>=3
                if bin==0 || mod(bin,1), error 'binning value must be a positive integer', end
                Z.bin = bin;
            end
            prepareFilter(Z,chgzoom,chgbin) % this will raise 'ChangedOperation' event
        end
        function moveZoom(Z,nstep)
            if strcmp(Z.zoom,':'), return, end
            d = diff(Z.zoom);
            if nstep > 0
                z2 = min(Z.headerin.n+.5,Z.zoom(2)+d*nstep);
                Z.setZoom([z2-d z2])
            else
                z1 = max(.5,Z.zoom(1)+d*nstep);
                Z.setZoom([z1 z1+d])
            end
        end
        function setBin(Z,bin)
            if bin==Z.bin, return, end
            % check
            if bin==0 || mod(bin,1), error 'binning value must be a positive integer', end
            % assign and update output
            Z.bin = bin;
            prepareFilter(Z,false,true) % this will raise 'ChangedOperation' event
        end
    end
    methods (Access='private')
        function prepareFilter(Z,chgzoom,chgbin)
            nin = Z.headerin.n;
            [curbin, curiout] = deal(Z.bin,Z.indicesout);
            
            % indices
            if strcmp(Z.zoom,':')
                if Z.bin==1
                    nout = nin;
                    Z.indicesin = 1:nin;
                else
                    nout = floor(nin/Z.bin);
                    Z.indicesin = reshape(1:Z.bin*nout,Z.bin,nout);
                end
            else
                if Z.bin==1
                    nout = fn_coerce(floor(1+diff(Z.zoom)), 1, nin);
                    idx1 = round(Z.zoom(1)+(1+diff(Z.zoom)-nout)/2);
                    Z.indicesin = idx1+(0:nout-1);
                else
                    nout = fn_coerce( floor(1+diff(Z.zoom)/Z.bin), 1, floor(nin/Z.bin) );
                    idx1 = round(Z.zoom(1)+(1+diff(Z.zoom)-nout*Z.bin)/2);
                    idx1 = fn_coerce(idx1,1,nin-Z.bin*nout+1);
                    Z.indicesin = reshape(idx1+(0:Z.bin*nout-1),Z.bin,nout);
                end
            end
            Z.indicesout = mean(Z.indicesin,1);
            
            % output header
            headin= Z.headerin;
            if strcmp(Z.zoom,':') && Z.bin==1
                % not any change
                Z.headerout = headin;
            elseif headin.ismeasure
                % if header is a measure, new positions are
                % straightforward to compute
                Z.headerout = xplr.header(headin.sublabels,nout,headin.start+(Z.indicesout(1)-1)*headin.scale,headin.scale*Z.bin);
            elseif headin.ncolumn==0
                % no values, keep track of index
                Z.headerout = xplr.header(headin.label,xplr.dimensionlabel('Index','numeric'),num2cell(Z.indicesout(:)));
            elseif Z.bin==1
                % no binning: getting values is straightforward
                Z.headerout = xplr.header(headin.label,headin.sublabels,headin.values(Z.indicesin,:));
            else
                % binning
                headvalue = Z.headerin.trackValues(num2cell(Z.indicesin,1));
                Z.headerout = xplr.header(headin.label,headin.sublabels,headvalue);
            end
            
            % notifications
            chgnout = (nout~=length(curiout));
            if chgzoom
                notify(Z,'ChangedOperation',xplr.eventinfo('zoom',chgnout))
            end
            if chgbin
                notify(Z,'ChangedOperation',xplr.eventinfo('bin'))
            end
            anychg = chgnout || (Z.indicesout(1)~=curiout(1));
            zoomin = anychg && (isempty(Z.indicesout) ...
                || ((Z.bin==1) && ~isempty(curiout) && (Z.indicesout(1)>=curiout(1)) && (Z.indicesout(end)<=curiout(end)) && (curbin==1)));
            if zoomin
                idxfirst = find(curiout==Z.indicesout(1),1,'first');
                idxlast = find(curiout==Z.indicesout(end),1,'last');
                idxrm = [1:idxfirst-1 idxlast+1:length(curiout)];
                notify(Z,'ChangedOperation',xplr.eventinfo('filter','remove',idxrm))
            elseif chgnout
                notify(Z,'ChangedOperation',xplr.eventinfo('filter','all'))
            elseif anychg
                notify(Z,'ChangedOperation',xplr.eventinfo('filter','chg',1:nout))
            end
        end
    end
    
    % Get Dependent
    methods
        function x = get.zoomvalue(Z)
            x = Z.zoom;
            if strcmp(x,':'), x = [.5 Z.headerin.n+.5]; end
        end
    end
    
    % Slicing
    methods
        function slic = slicing(Z,dat,dims,selsubidx)
            % here Z can be non-scalar!
            if length(dims)~=length(Z), error 'number of dimensions does not match number of points', end
            dosubidx = (nargin>=4);
            
            % size
            s = size(dat);
            nddata = max(max(dims),length(s));
            s(end+1:nddata) = 1;
            sout = s;
            headout = [Z.headerout];
            sout(dims) = [headout.n];
            
            % extract sub-data
            subs = substruct('()',repmat({':'},1,length(s)));
            nofilt = true;
            for i=1:length(Z)
                if ~strcmp(Z(i).zoom,':') || Z(i).bin~=1 || dosubidx
                    nofilt = false;
                    ind = Z(i).indicesin;
                    if dosubidx, ind = ind(:,selsubidx); end
                    subs.subs{dims(i)} = ind; 
                end
            end
            if nofilt
                slic = dat;
            else
                slic = subsref(dat,subs);
            end
            
            % bin
            scur = size(slic); scur(end+1:nddata) = 1; % size before binning
            for i=1:length(Z)
                if Z(i).bin>1
                    slic = reshape(slic,[prod(sout(1:dims(i)-1)) Z(i).bin sout(dims(i)) prod(scur(dims(i)+1:end))]);
                    slic = nmean(slic,2);
                end
            end
            slic = reshape(slic,sout);
        end
    end
    methods (Access='protected')
        function slic = operation_(Z,dat,dims)
            % function slic = operation_(Z,dat,dims)
            %---
            % dat and slic are simple Matlab arrays
            slic = Z.slicing(dat,dims);
        end
        function updateOperation_(Z,x,dims,slice,flag,ind) 
            % there is no 'smart' way of updating operation: just do the
            % slicing
            slic = Z.slicing(x.data,dims);
            slice.updateData(flag,dims,ind,slic,Z.headerout); % this will trigger automatic notifications
        end
    end
    
    % Link with zoom definition in real world coordinates
    methods
        function zoom_world = operationData2Space(Z)
            if strcmp(Z.zoom,':')
                zoom_world = ':';
            else
                zoom_world = Z.headerin.start + (Z.zoom-1)*Z.headerin.scale;
            end
        end
        function updateOperationData2Space(Z,WO,evnt)
            if ~strcmp(evnt.type,'zoom'), return, end
            WO.operation = Z.operationData2Space();
            notify(WO,'ChangedOperation')
        end
        function updateOperationSpace2Data(Z,world_operation,~)
            if strcmp(world_operation,':')
                Z.setZoom(':')
            else
                zoom_idx = 1 + (world_operation - Z.headerin.start)/Z.headerin.scale;
                Z.setZoom(zoom_idx)
            end
        end
    end
    
    % Tools
    methods
        function idx1 = orig2zoomed(Z,idx)
            b = (idx>=Z.indicesin(1)) & (idx<=Z.indicesin(end));
            idx1 = idx;
            idx1(~b) = 0;
            idx1(b) = 1 + floor((idx(b)-Z.indicesin(1))/Z.bin);
        end
    end
    
end