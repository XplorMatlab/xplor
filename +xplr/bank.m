classdef bank < hgsetget
    
    % Content
    properties (SetAccess='private')
        % views
        currentviews = struct('obj',cell(1,0),'hl',cell(1,0)); % all current views are registered here
        % units
        measures = struct('label','time','units',struct('unit',{'s' 'ms'},'value',{1 1e-3}));
        % headers
        recentheaders = xplr.header.empty(1,0); % headers will be ordered according to their appearance date
        % filter sets
        filtersets = xplr.filterSet.empty(1,0);
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
            xplr.bank.registerheaders(V.data.header)
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
            if exist(fsave,'file')
                save(fsave,'-STRUCT','s','-APPEND');
            else
                save(fsave,'-STRUCT','s');
            end
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
    methods (Static)
        function registerheaders(newheader)
            B = xplr.bank.getbank();
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
        function head = getrecentheaders(n, num_max)
            B = xplr.bank.getbank();
            % get a list of recent headers for data length n
            same_length = [B.recentheaders.n]==n;
            head = B.recentheaders(same_length);
            % add also recent enumeration headers, even though they might
            % be of different length
            enum = B.recentheaders(~same_length & [B.recentheaders.isenum]);
            [~, idx] = unique({enum.label},'stable');
            head = [head enum(idx)];
            if nargin>=2, head(num_max+1:end) = []; end
        end
        function clearrecentheaders()
            B = xplr.bank.getbank();
            B.recentheaders(:) = [];
            saveprop(B,'recentheaders')
        end
    end
    
    % Filter sets
    methods (Static, Access=private)
        function FS = getFilterSet(linkkey)
            B = xplr.bank.getbank();
            if linkkey<=length(B.filtersets)
                FS = B.filtersets(linkkey);
            elseif linkkey==length(B.filtersets)+1
                FS = xplr.filterSet(linkkey);
                B.filtersets(end+1) = FS;
            else
                error('New filter set key should be %i, encountered %i instead.',length(B.filtersets)+1,idx)
            end
        end
    end
    methods (Static)
        function n = nfilterset()
            B = xplr.bank.getbank();
            n = length(B.filtersets);
        end
        function F = getFilter(linkkey, head, varargin)
            % function F = getFilter(linkkey, header[,user])
            FS = xplr.bank.getFilterSet(linkkey);
            F = FS.getFilter(head, varargin{:});
        end
        function addFilter(linkkey, F, varargin)
            % function addFilter(linkkey, F[,user])
            FS = xplr.bank.getFilterSet(linkkey);
            FS.addFilter(F, varargin{:})
        end
        function removeFilter(linkkey, F, user)
            % function removeFilter(linkkey, F, user)
            FS = xplr.bank.getFilterSet(linkkey);
            FS.removeFilter(F, user)
        end            
        function showList(linkkey, F, user)
            % function removeFilter(linkkey, F, user)
            FS = xplr.bank.getFilterSet(linkkey);
            FS.removeFilter(F, user)
        end
        function F = getZoomFilter(linkkey, head, varargin)
            % function F = getZoomFilter(linkkey, header[,user])
            FS = xplr.bank.getFilterSet(linkkey);
            F = FS.getZoomFilter(head, varargin{:});
        end
        function addZoomFilter(linkkey, F, varargin)
            % function addZoomFilter(linkkey, F[,user])
            for i=F
                FS = xplr.bank.getFilterSet(linkkey);
                FS.addZoomFilter(i, varargin{:})
            end
        end
        function removeZoomFilter(linkkey, F, user)
            % function removeZoomFilter(linkkey, F, user)
            FS = xplr.bank.getFilterSet(linkkey);
            FS.removeZoomFilter(F, user)
        end            
    end
end