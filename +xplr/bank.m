classdef bank < handle
% xplr.bank The bank is unique and store several elements :
%
%  * currentviews: not used yet but it will be used for multiview
%  * measures: units for measures conversion
%  * recent headers: previously used headers, this list of recent headers is scanned when a new set of data is xplored
%  * filters_registry: registry of all existing filters, indexed by filter type, link key and input header
%  * list_combo: xplr.listcombo object for displaying 1D filters as lists

    properties (SetAccess='private')
        currentviews = struct('obj',cell(1,0),'hl',cell(1,0));                                  % all current views are registered here
        measures = struct('label','time','units',struct('unit',{'s' 'ms'},'value',{1 1e-3}));   % units
        recentheaders = xplr.header.empty(1,0);                                                 % headers will be ordered according to their appearance date         
        filters_registry
        list_combo
    end
    
    % Constructor
    methods (Access='private')
        function B = bank()
            % bank constructor
            
            % filters filters_registry
            B.filters_registry = xplr.bankRegistry;
            
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
        % This method is static because there will be one specialized
        % method by filter type.
        function F = getExistingFilter(filtertype, linkkey, header, user)
            % function F = getRegistryValue(filtertype, linkkey, header[, newuser])
            B = xplr.bank.getbank();
            hID = getID(header);
            F = B.filters_registry.getValue({filtertype, linkkey, hID},user);
            % F was in fact deleted? -> unregister
            if ~isempty(F) && ~isvalid(F)
                B.filters_registry.unregister({filtertype, linkkey, hID});
            end
        end
        function [F, isnew] = getFilter(filtertype, linkkey, header, user)
            if strcmp(filtertype,'filterAndPoint')
                error 'getting a filterAndPoint filter necessitates a specialized method, so calling getFilter is not authorized'
            end
            F = xplr.bank.getExistingFilter(filtertype, linkkey, header, user);
            if isempty(F)
                isnew = true;
                F = feval(['xplr.' filtertype],header);
                xplr.bank.registerFilter(linkkey,F,user);
                % if input space is measurable, connect with a
                % 'worldOperand' object that will synchronize several filters
                % corresponding to different referentials in this space
                key = header.getMeasureSpaceID;
                if ~isempty(key)
                   % search a corresponding worldOperand in the filters
                   % registry, create if it does not exist
                   B = xplr.bank.getbank();
                   worldfiltertype = ['world' filtertype];
                   worldfilter = B.filters_registry.getValue({worldfiltertype, linkkey, key}, F);
                   if isempty(worldfilter)
                       % create new worldOperand object to link to F
                       worldfilter = xplr.worldOperand(F);
                       B.filters_registry.register({worldfiltertype, linkkey, key}, worldfilter, F);
                   else
                       worldfilter.connectDataOperand(F);
                   end
                end
            else
                isnew = false;
            end
        end
    end
    methods (Static)
        function keys = availableFilterKeys(filtertype)
            B = xplr.bank.getbank();
            % We recall the nested structure of filters registry: filters
            % are indexed first by filter type, then by link key, then by
            % header ID. 
            % First get the sub-registry corresponding to filter type.
            sub_registry = B.filters_registry.getValue(filtertype);
            % This sub-registry is indexed by link keys.
            if isempty(sub_registry)
                keys = zeros(1,0);
            else
                keys = [sub_registry.content.key];
            end
        end
        function registerFilter(linkkey, F, user)
            % function registerFilter(linkkey, F, user)
            B = xplr.bank.getbank();
            for Fi = row(F)
                filtertype = strrep(class(Fi),'xplr.','');
                hID = getID(Fi.headerin);
                Fi.linkkey = linkkey;   % memorize linkkey inside filter
                B.filters_registry.register({filtertype, linkkey, hID}, Fi, user)
            end
        end
        function unregisterFilter(F, user)
            % function unregisterFilter(F, user)
            B = xplr.bank.getbank();
            for Fi = row(F)
                filtertype = strrep(class(Fi),'xplr.','');
                if ~isvalid(Fi), continue, end
                linkkey = Fi.linkkey;
                hID = getID(Fi.headerin);
                B.filters_registry.unregister({filtertype, linkkey, hID}, user);
            end
        end
        % Below are the specialized methods for getting filters,
        % specialized by filter type. They create the filter if it does not
        % exist. This creation sometimes relies on accessing/creating
        % filters of another type (for example filterAndPoint filters
        % use/create a filter and a point filter; zoomfilter filters
        % use/create a zoomcentral object)
        function F = getFilterAndPoint(linkkey, header, user, doshow)
            % function F = getFilterAndPoint(linkkey, header, newuser[, doshow])
            % function F = getFilterAndPoint([linkkey_filter linkkey_point], header[, doshow[, newuser]])
            if nargin<4, doshow = false; end
            F = xplr.bank.getExistingFilter('filterAndPoint', linkkey, header, user);
            if isempty(F)
                % construct filterAndPoint object from filter and point
                % objects obtained themselves from the bank
                FF = xplr.bank.getFilterFilter(linkkey, header, user);
                FP = xplr.bank.getPointFilter(linkkey, header, user);
                F = xplr.filterAndPoint(FF,FP);
                xplr.bank.registerFilter(linkkey,F,user);
                if doshow 
                    if F.ndin > 1
                        disp 'cannot display list for ND filter'
                    else
                        xplr.bank.showList(F)
                    end
                end     
            end
        end
        function F = getFilterFilter(linkkey, header, user)
            % function F = getFilterFilter(linkkey, header[, newuser]])
            F = xplr.bank.getFilter('filter', linkkey, header, user);
        end
        function P = getPointFilter(linkkey, header, user)
            % function F = getPointFilter(linkkey, header[, newuser]])
            for i = 1:length(header)
                P(i) = xplr.bank.getFilter('point', linkkey, header(i), user);
            end
        end
        function F = getZoomFilter(linkkey, header, user)
            % function F = getZoomFilter(B,header, user)
            F = xplr.bank.getFilter('zoomfilter', linkkey, header, user);
        end
        function showList(F)
            if ~isa(F,'xplr.filterAndPoint')
                error 'only filterAndPoint object can be shown'
            elseif F.ndin > 1
                % not possible to show list if filter is not 1D, ignore
                return
            end
            
            B = xplr.bank.getbank();
            % Create list combo?
            if isempty(B.list_combo) || ~isvalid(B.list_combo)
                B.list_combo = xplr.listcombo();
%                 % no need to delete the listener upon filterSet deletion: filterSet are supposed never to be deleted
%                 connectlistener(B.list_combo,B,'Empty',@(u,e)set(B,'list_combo',[])); 
            end
            
            B.list_combo.showList(F)
            figure(B.list_combo.container.hobj)
        end
    end
end