classdef Bank < handle
% xplr.Bank The bank is unique and store several elements :
%
%  * measures: units for measures conversion
%  * recent headers: previously used headers, this list of recent headers is scanned when a new set of data is xplored
%  * filters_registry: registry of all existing filters, indexed by filter type, link key and input header
%  * list_combo: xplr.listcombo object for displaying 1D filters as lists

    properties (SetAccess='private')
        measures = struct('label', 'time', 'units', struct('unit', {'s' 'ms'}, 'value', {1 1e-3}));   % units
        recent_headers = xplr.Header.empty(1,0);                                                 % headers will be ordered according to their appearance date
        filters_registry
        list_combo
    end
    
    % Constructor
    methods (Access='private')
        function B = Bank()
            % bank constructor
            
            % filters filters_registry
            B.filters_registry = xplr.BankRegistry;
            
            % load saved measures
            B.load_prop('measures')
            
            % load saved recent headers
            B.load_prop('recent_headers')
            B.recent_headers(~isvalid(B.recent_headers)) = [];
        end
    end
    methods (Static)
        function B = get_bank()
            % Unique bank object is stored in a persistent variable.
            persistent B0
            if isempty(B0)
                B0 = xplr.Bank();
            end
            B = B0;
        end
    end
    
    % Views
    methods (Static)
        function register_view(V)
            % update list of recent headers
            dim_header = V.data.header; % xplr.dim_header class, need to convert to xplr.Header!
            xplr.Bank.register_headers(xplr.Header(dim_header))
        end
    end
    
    % Load/save field
    methods (Access='private')
        function load_prop(B, prop)
            % load_prop
            f_save = fn_userconfig('configfolder', 'xplor_bank');
            warning('off', 'MATLAB:load:variableNotFound')
            try %#ok<TRYNC>
                B.(prop) = fn_loadvar(f_save, prop);
            end
            warning('on','MATLAB:load:variableNotFound')
        end
        function save_prop(B, prop)
            f_save = fn_userconfig('configfolder', 'xplor_bank');
            s = struct(prop, {B.(prop)});
            if exist(f_save, 'file')
                save(f_save, '-STRUCT', 's', '-APPEND');
            else
                save(f_save, '-STRUCT', 's');
            end
        end
    end
    
    % Units
    methods (Static)
        function [measure_label, conversion, measure] = get_unit_info(unit)
            % function [measurelabel, conversion, measure] = get_unit_info(unit)
            %
            % if 'unit' is a registered unit (e.g. 'ms'), returns the label 
            % of the corresponding measure (e.g. 'time') and the conversion
            % to the reference unit (e.g. 1e-3, reference unit being 's') 
            % otherwise returns an empty label
            B = xplr.Bank.get_bank();
            for i=1:length(B.measures)
                mi = B.measures(i);
                idx = strcmp({mi.units.unit}, unit);
                if any(idx)
                    measure_label = mi.label;
                    conversion = mi.units(idx).value;
                    measure = mi;
                    return
                end
            end
            [measure_label, conversion, measure] = deal([]);
        end
        function m = get_measures()
            m = xplr.Bank.get_bank().measures;
        end
        function add_measure(label, units)
            % check
            if ~ischar(label) || ~isstruct(units) || ~isequal(fieldnames(units), {'unit'; 'value'})
                error 'new measure is not properly defined'
            end
            % set
            B = xplr.Bank.get_bank();
            if any(strcmp({B.measures.label}, label))
                error('a measure labelled ''%s'' already exists', label)
            end
            B.measures(end+1) = struct('label',label,'units',units);
            B.save_prop('measures')
        end
        function edit_measure(oldlabel, label, units)
            % check
            if ~ischar(label) || ~isstruct(units) || ~isequal(fieldnames(units),{'unit'; 'value'})
                error 'new measure is not properly defined'
            end
            % set
            B = xplr.Bank.get_bank();
            idx = strcmp({B.measures.label}, oldlabel);
            if ~any(idx), error('no existing measure is labelled ''%s''', label), end
            B.measures(idx) = struct('label', label, 'units', units);
            B.save_prop('measures')
        end
    end
    
    % Headers
    methods (Static)
        function register_headers(new_header)
            B = xplr.Bank.get_bank();
            new_header(fn_isemptyc({new_header.label})) = [];
            n = length(new_header);
            idx = cell(1, n);
            for i = 1:n
                idx{i} = fn_find(new_header(i), B.recent_headers, 'first');
            end
            idx = [idx{:}];
            if isequal(idx, 1:n), return, end % new headers are already at the beginning of the list
            % place all new headers first in the list
            B.recent_headers(idx(idx ~= 0)) = [];
            B.recent_headers = [new_header, B.recent_headers];
            n_header_max = xplr.Parameters.get('bank.n_header_max');
            B.recent_headers(n_header_max+1:end) = [];
            save_prop(B,'recent_headers')
        end
        function head = get_recent_headers(n, num_max)
            B = xplr.Bank.get_bank();
            % get a list of recent headers for data length n
            same_length = [B.recent_headers.n] == n;
            head = B.recent_headers(same_length);
            % add also recent enumeration headers, even though they might
            % be of different length
            enum = B.recent_headers(~same_length & [B.recent_headers.is_enum]);
            [~, idx] = unique({enum.label}, 'stable');
            head = [head, enum(idx)];
            if nargin >= 2, head(num_max+1:end) = []; end
        end
        function clear_recent_headers()
            B = xplr.Bank.get_bank();
            B.recent_headers(:) = [];
            save_prop(B, 'recent_headers')
        end
    end
    
    % Filter sets
    methods (Static, Access='private')
        % This method is static because there will be one specialized
        % method by filter type.
        function F = get_existing_filter(filter_type, link_key, header, user)
            % function F = get_existing_filter(filter_type, link_key, header[, newuser])
            B = xplr.Bank.get_bank();
            header = xplr.Header(header); % in case header is a xplr.dimheader
            h_id = get_id(header);
            F = B.filters_registry.get_value({filter_type, link_key, h_id}, user);
            % F was in fact deleted? -> unregister
            if ~isempty(F) && ~isvalid(F)
                B.filters_registry.unregister({filter_type, link_key, h_id});
            end
        end
        function [F, is_new] = get_filter(filter_type, link_key, header, user)
            if strcmp(filter_type,'FilterAndPoint')
                error 'getting a FilterAndPoint filter necessitates a specialized method, so calling getFilter is not authorized'
            end
            header = xplr.Header(header); % in case header is a xplr.dimheader
            F = xplr.Bank.get_existing_filter(filter_type, link_key, header, user);
            if isempty(F)
                isnew = true;
                F = feval(['xplr.', filter_type], header);
                xplr.Bank.register_filter(link_key, F, user);
                % if input space is measurable, connect with a
                % 'worldOperand' object that will synchronize several filters
                % corresponding to different referentials in this space
                key = header.get_measure_space_id;
                if ~isempty(key)
                   % search a corresponding worldOperand in the filters
                   % registry, create if it does not exist
                   B = xplr.Bank.get_bank();
                   world_filter_type = ['world', filter_type];
                   world_filter = B.filters_registry.get_value({world_filter_type, link_key, key}, F);
                   if isempty(world_filter)
                       % create new worldOperand object to link to F
                       world_filter = xplr.WorldOperand(F);
                       B.filters_registry.register({world_filter_type, link_key, key}, world_filter, F);
                   else
                       world_filter.connect_data_operand(F);
                   end
                end
            else
                is_new = false;
            end
        end
    end
    methods (Static)
        function keys = available_filter_keys(filter_type)
            B = xplr.Bank.get_bank();
            % We recall the nested structure of filters registry: filters
            % are indexed first by filter type, then by link key, then by
            % header ID. 
            % First get the sub-registry corresponding to filter type.
            sub_registry = B.filters_registry.get_value(filter_type);
            % This sub-registry is indexed by link keys.
            if isempty(sub_registry)
                keys = zeros(1, 0);
            else
                keys = [sub_registry.content.key];
            end
        end
        function register_filter(link_key, F, user)
            % function register_filter(link_key, F, user)
            B = xplr.Bank.get_bank();
            for Fi = row(F)
                filter_type = strrep(class(Fi), 'xplr.', '');
                h_id = get_id(Fi.header_in);
                Fi.link_key = link_key;   % memorize linkkey inside filter
                B.filters_registry.register({filter_type, link_key, h_id}, Fi, user)
            end
        end
        function unregister_filter(F, user)
            % function unregister_filter(F, user)
            B = xplr.Bank.get_bank();
            for Fi = row(F)
                filter_type = strrep(class(Fi), 'xplr.', '');
                if ~isvalid(Fi), continue, end
                link_key = Fi.link_key;
                h_id = get_id(Fi.header_in);
                B.filters_registry.unregister({filter_type, link_key, h_id}, user);
            end
        end
        % Below are the specialized methods for getting filters,
        % specialized by filter type. They create the filter if it does not
        % exist. This creation sometimes relies on accessing/creating
        % filters of another type (for example filterAndPoint filters
        % use/create a filter and a point filter; zoomfilter filters
        % use/create a zoomcentral object)
        function F = get_filter_and_point(link_key, header, user, do_show)
            % function F = get_filter_and_point(link_key, header, newuser[, do_show])
            % function F = get_filter_and_point([link_key_filter linkkey_point], header[, do_show[, new_user]])
            if nargin<4, do_show = false; end
            F = xplr.Bank.get_existing_filter('FilterAndPoint', link_key, header, user);
            if isempty(F)
                % construct filterAndPoint object from filter and point
                % objects obtained themselves from the bank
                % (get filter and point from the bank, do not register a
                % user yet)
                FF = xplr.Bank.get_filter_filter(link_key, header, []);
                FP = xplr.Bank.get_point_filter(link_key, header, []);
                % (create FilterAndPoint object)
                F = xplr.FilterAndPoint(FF,FP);
                % (now register F as user of FF and FP)
                xplr.Bank.get_filter_filter(link_key, header, F);
                xplr.Bank.get_point_filter(link_key, header, F);
                % (and register user as a user of F)
                xplr.Bank.register_filter(link_key, F, user);
                % show filter
                if do_show
                    if F.nd_in > 1
                        disp 'cannot display list for ND filter'
                    else
                        xplr.Bank.show_list(F)
                    end
                end     
            end
        end
        function F = get_filter_filter(link_key, header, user)
            % function F = get_filter_filter(link_key, header[, new_user]])
            F = xplr.Bank.get_filter('Filter', link_key, header, user);
        end
        function P = get_point_filter(link_key, header, user)
            % function F = get_point_filter(link_key, header[, new_user]])
            for i = 1:length(header)
                P(i) = xplr.Bank.get_filter('Point', link_key, header(i), user);
            end
        end
        function F = get_zoom_filter(link_key, header, user)
            % function F = get_zoom_filter(B, header, user)
            F = xplr.Bank.get_filter('ZoomFilter', link_key, header, user);
        end
        function show_list(F)
            if ~isa(F,'xplr.FilterAndPoint')
                error 'only FilterAndPoint object can be shown'
            elseif F.nd_in > 1
                % not possible to show list if filter is not 1D, ignore
                return
            end
            
            B = xplr.Bank.get_bank();
            % Create list combo?
            if isempty(B.list_combo) || ~isvalid(B.list_combo)
                B.list_combo = xplr.ListCombo();
%                 % no need to delete the listener upon filterSet deletion: filterSet are supposed never to be deleted
%                 connect_listener(B.list_combo, B, 'Empty', @(u,e)set(B, 'list_combo', []));
            end
            
            B.list_combo.show_list(F)
            figure(B.list_combo.container.hobj)
        end
    end
end