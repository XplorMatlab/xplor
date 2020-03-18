classdef bank < hgsetget
% xplr.bank The bank is unique and store several elements :
%
%  * currentviews: not used yet but it will be used for multiview
%  * recent headers: previously used headers, this list of recent headers is scanned when a new set of data is xplored
%  * measures: units for measures conversion
%  * filterssets: filters and zoomfilters stored by linkkey
%   * linkey
%   * registry
%   * combo
%   * zregistry
%   * zcregistry
%

    properties (SetAccess='private')
        currentviews = struct('obj',cell(1,0),'hl',cell(1,0));                                  % all current views are registered here
        measures = struct('label','time','units',struct('unit',{'s' 'ms'},'value',{1 1e-3}));   % units
        recentheaders = xplr.header.empty(1,0);                                                 % headers will be ordered according to their appearance date 
        filtersets = xplr.filterSet.empty(1,0);                                                 % filter sets
    end
    
    methods (Access='private')
        function B = bank()
            % bank constructor

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
            dimheader = V.data.header; % xplr.dimheader class, need to convert to xplr.header!
            xplr.bank.registerheaders(xplr.header(dimheader))
        end
    end
    
    % Load/save field
    methods (Access='private')
        function loadprop(B,prop)
            % loadprop
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
            %
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
        function unregisterView(V)
            % unregister view
            B = xplr.bank.getbank();
            idx = ([B.currentviews.obj]==V);
            delete([B.currentviews(idx).hl])
            B.currentviews(idx) = [];
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
    methods (Static, Access='private')
        function FS = getFilterSet(linkkey)
            B = xplr.bank.getbank();
            if linkkey<=length(B.filtersets)
                FS = B.filtersets(linkkey);
            elseif linkkey==length(B.filtersets)+1
                FS = xplr.filterSet(linkkey);
                B.filtersets(end+1) = FS;
            else
                error('New filter set key should be %i, encountered %i instead.',length(B.filtersets)+1,linkkey)
            end
        end
    end
    methods (Static)
        function n = nfilterset()
            B = xplr.bank.getbank();
            n = length(B.filtersets);
        end
        function F = getFilter(linkkey, head, doshow, varargin)
            % function F = getFilter(linkkey, header ,doshow [,user])
            if nargin<3, doshow = false; end
            FS = xplr.bank.getFilterSet(linkkey);
            F = FS.getFilter(head, doshow, varargin{:});
        end
        function registerFilter(linkkey, F, varargin)
            % function registerFilter(linkkey, F[,user])
            FS = xplr.bank.getFilterSet(linkkey);
            doshow=true;
            FS.addFilter(F, doshow, varargin{:})
        end
        function unregisterFilter(F, user)
            % function unregisterFilter(F, user)
            if F.linkkey == 0
                disp 'attempt to unregister a private filter from the bank!'
                return
            end            
            FS = xplr.bank.getFilterSet(F.linkkey);
            FS.removeFilter(F, user)
        end            
        function showList(F)
            % function showList(F)
            if F.linkkey == 0
                disp 'attempt to show the list display from the bank for a private filter!'
                return
            end            
            FS = xplr.bank.getFilterSet(F.linkkey);
            FS.showList(F)
        end
        function F = getZoomFilter(linkkey, head, varargin)
            % function F = getZoomFilter(linkkey, header[,user])
            FS = xplr.bank.getFilterSet(linkkey);
            F = FS.getZoomFilter(head, varargin{:});
        end
        function registerZoomFilter(linkkey, F, varargin)
            % function registerZoomFilter(linkkey, F[,user])
            for i=F
                FS = xplr.bank.getFilterSet(linkkey);
                FS.addZoomFilter(i, varargin{:})
            end
        end
        function unregisterZoomFilter(F, user)
            % function unregisterZoomFilter(F, user)
            if F.linkkey == 0
                disp 'attempt to unregister a private zoom filter from the bank!'
                return
            end            
            FS = xplr.bank.getFilterSet(F.linkkey);
            FS.removeZoomFilter(F, user)
        end       
        function keys = availableFilterKeys()
            B = xplr.bank.getbank();
            keys = [B.filtersets.linkkey];
        end
    end
end