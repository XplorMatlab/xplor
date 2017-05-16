
classdef bank < hgsetget
    
    % Content
    properties (SetAccess='private')
        % views
        currentviews = struct('obj',cell(1,0),'hl',cell(1,0)); % all current views are registered here
        % units
        measures = struct('label','time','units',struct('unit',{'s' 'ms'},'value',{1 1e-3}));
        % headers
        recentheaders = xplr.header.empty(1,0); % headers will be ordered according to their appearance date
        %         filters = struct('hID',cell(1,0),'key',[],'obj',[],'indicesonly',[]);
        %         points = struct('hID',cell(1,0),'key',[],'obj',[]);
    end
    
    % (old)
    methods (Access='private')
        %         function [keyavail idxavail keyforbid] = checkkeys(B,type,hID,indicesonly)
        %             switch type
        %                 case 'point'
        %                     % no restriction exist
        %                     idxavail = find([B.points.hID]==sum(hID));
        %                     keyavail = [B.points(idxavail).key];
        %                     keyforbid = [];
        %                 case 'filter'
        %                     % restriction: there must be no overlap in header IDs
        %                     % with other existing filters
        %                     idxavail = find(fn_map(@sum,{B.filters.hID})+[B.filters.indicesonly]==sum(hID)+indicesonly);
        %                     keyavail = [B.filters(idxavail).key];
        %                     idxintersect = fn_find({B.filters.hID},@(id)any(ismember(id,hID)),'all');
        %                     idxforbid = setdiff(idxintersect,idxavail);
        %                     keyforbid = unique([B.filters(idxforbid).key]);
        %             end
        %         end
        %         function [obj key] = getobject(B,type,headerin,key,indicesonly)
        %             if nargin<5, indicesonly = false; end
        %             hID = [headerin.ID];
        %             [keyavail idxavail keyforbid] = checkkeys(B,type,hID,indicesonly);
        %             if ~isempty(key) && ismember(key,keyforbid)
        %                 error('no filter is available for these headers and with key %i',key)
        %             end
        %             if isempty(keyavail)
        %                 key = 0;
        %                 while ismember(key,keyforbid), key = key+1; end
        %                 switch type
        %                     case 'point'
        %                         obj = xplr.xpoint(headerin);
        %                         B.points(end+1) = struct('hID',hID,'key',key,'obj',obj);
        %                     case 'filter'
        %                         obj = xplr.xfilter(headerin,indicesonly);
        %                         B.filters(end+1) = struct('hID',hID,'key',key,'obj',obj,'indicesonly',indicesonly);
        %                 end
        %             else
        %                 [~, j] = min(keyavail);
        %                 idx = idxavail(j);
        %                 switch type
        %                     case 'point'
        %                         obj = B.points(idx).obj;
        %                     case 'filter'
        %                         obj = B.filters(idx).obj;
        %                 end
        %             end
        %         end
    end
    
    % Bank
    methods (Access='private')
        function B = bank()
            % load saved measures
            B.loadprop('measures')
            % load saved recent headers
            B.loadprop('recentheaders')
        end
    end
    methods (Static)
        function B = getbank()
            % Unique filter bank is attached to the root graphic object.
            % This is preferrable to using a global variable that might be
            % deleted with the 'clear' command.
%             if isappdata(0,'xplr_filterbank')
%                 B = getappdata(0,'xplr_filterbank');
%             else
%                 B = xplr.bank();
%                 setappdata(0,'xplr_filterbank',B)
%             end
            persistent B0
            if isempty(B0)
                B0 = xplr.bank();
            end
            B = B0;
        end
    end
    
    % Views
    methods (Static)
        function registerView(V)
            % register view
            B = xplr.bank.getbank();
            hl = addlistener(V,'ObjectBeingDestroyed',@(u,e)xplr.bank.unregisterView(V));
            B.currentviews(end+1) = struct('obj',V,'hl',hl);
            % update list of recent headers
            B.registerheadersPrivate(V.data.header)
        end
        function unregisterView(V)
            B = xplr.bank.getbank();
            idx = ([B.currentviews.obj]==V);
            delete([B.currentviews(idx).hl])
            B.currentviews(idx) = [];
        end
    end
    
    % Load/save field
    methods (Access='private')
        function loadprop(B,prop)
            fsave = fn_userconfig('configfolder','xplr.bank');
            warning('off','MATLAB:load:variableNotFound')
            try %#ok<TRYNC>
                B.(prop) = fn_loadvar(fsave,prop);
            end
            warning('on','MATLAB:load:variableNotFound')
        end
        function saveprop(B,prop)
            fsave = fn_userconfig('configfolder','xplr.bank');
            s = struct(prop,B.(prop)); %#ok<NASGU>
            save(fsave,'-STRUCT','s','-APPEND');
        end
    end
    
    % Units
    methods (Static)
        function [measurelabel conversion measure] = getunitinfo(unit)
            % function [measurelabel conversion measure] = getunitinfo(unit)
            %---
            % if 'unit' is a registered unit (e.g. 'ms'), returns the label 
            % of the corresponding measure (e.g. 'time') and the conversion
            % to the reference unit (e.g. 1e-3, reference unit being 's') 
            % otherwise returns an empty label
            B = xplr.bank.getbank();
            for i=1:length(B.measures)
                mi = B.measures(i);
                idx = strcmp({mi.units.unit},unit);
                if any(idx)
                    measurelabel = mi.label;
                    conversion = mi.units(idx).value;
                    measure = mi;
                    return
                end
            end
            [measurelabel conversion measure] = deal([]);
        end
        function m = getMeasures()
            m = xplr.bank.getbank().measures;
        end
        function addMeasure(label,units)
            % check
            if ~ischar(label) || ~isstruct(units) || ~isequal(fieldnames(units),{'unit'; 'value'})
                error 'new measure is not properly defined'
            end
            % set
            B = xplr.bank.getbank();
            if any(strcmp({B.measures.label},label))
                error('a measure labelled ''%s'' already exists',label)
            end
            B.measures(end+1) = struct('label',label,'units',units);
            B.saveprop('measures')
        end
        function editMeasure(oldlabel,label,units)
            % check
            if ~ischar(label) || ~isstruct(units) || ~isequal(fieldnames(units),{'unit'; 'value'})
                error 'new measure is not properly defined'
            end
            % set
            B = xplr.bank.getbank();
            idx = strcmp({B.measures.label},oldlabel);
            if ~any(idx), error('no existing measure is labelled ''%s''',label), end
            B.measures(idx) = struct('label',label,'units',units);
            B.saveprop('measures')
        end
    end
    
    % Headers
    methods (Access='private')
        function registerheadersPrivate(B,newheader)
            newheader(fn_isemptyc({newheader.label})) = [];
            n = length(newheader);
            idx = cell(1,n);
            for i=1:n
                idx{i} = fn_find(newheader(i),B.recentheaders,'first'); 
            end
            idx = [idx{:}];
            if isequal(idx,1:n), return, end % new headers are already at the beginning of the list
            % place all new headers first in the list
            B.recentheaders(idx(idx~=0)) = [];
            B.recentheaders = [newheader B.recentheaders];
            nheadermax = xplr.parameters.get('bank.nheadermax');
            B.recentheaders(nheadermax+1:end) = [];
            saveprop(B,'recentheaders')
        end
    end
    methods (Static)
        function registerheaders(newheader)
            xplr.bank.getbank.registerheadersPrivate(newheader)
        end
        function head = getrecentheaders(n)
            B = xplr.bank.getbank();
            % get a list of recent headers for data length n
            head = B.recentheaders([B.recentheaders.n]==n);
        end
        function clearrecentheaders()
            B = xplr.bank.getbank();
            B.recentheaders(:) = [];
            saveprop(B,'recentheaders')
        end
    end
    
    % Filters (old)
    methods (Static)
        %         function [F key] = getfilter(headerin,key,indicesonly)
        %             if ~isscalar(headerin), error 'header input must be scalar', end
        %             if nargin<2, key = []; end
        %             if nargin<3, indicesonly = false; end
        %             B = xplr.bank.getbank();
        %             [F key] = B.getobject('filter',headerin,key,indicesonly);
        %         end
        %         function [P key] = getpoint(headerin,key)
        %             if nargin<2, key = []; end
        %             B = xplr.bank.getbank();
        %             for i=1:length(headerin)
        %                 P(i) = B.getobject('point',headerin,key); %#ok<AGROW>
        %             end
        %         end
        %         function [F P key] = getFilterAndPoint(headerin,key,indicesonly)
        %             if ~isscalar(headerin), error 'header input must be scalar', end
        %             if nargin<2, key = []; end
        %             if nargin<3, indicesonly = false; end
        %             [F key] = xplr.bank.getfilter(headerin,key,indicesonly);
        %             P = xplr.bank.getpoint(headerin,key);
        %         end
    end
    
end