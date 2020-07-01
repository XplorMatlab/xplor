classdef DataOperand < xplr.GraphNode
% dataOperand
% Abstract class defining an operation on an xdata object
    
    properties (SetAccess='protected')
        header_in
        header_out
    end
    properties (Dependent, SetAccess='private', Transient)
        sz_in
        nd_in
        sz_out
        nd_out
        reduction_factor
    end
    % properties below are not handled by dataOperand class and
    % sub-classes, but rather by the objects that use them; they should not
    % be set by user however!!
    properties (Transient)
        link_key = 0
        world_operand
    end
   
    % There are two events: operation definition can change without the
    % operation itself being changed, for example when only slightly moving
    % a time cursor its position (operation definition) has changed, but
    % not the pixel it selects (underlying slicing operation is unchanged).
    events
        ChangedOperation
    end
    
    % Constructor
    methods
        
    end
    
    % Size
    methods
        function x = get.reduction_factor(O)
            x = prod([O.header_in.n])/prod([O.header_out.n]);
        end
        function sz = get.sz_in(O)
            sz = [O.header_in.n];
        end
        function nd = get.nd_in(O)
            nd = length(O.header_in);
        end
        function sz = get.sz_out(O)
            sz = [O.header_out.n];
        end
        function nd = get.nd_out(O)
            nd = length(O.header_out);
        end
    end
    
    % Operation
    methods (Abstract, Access='protected')
        dat = operation_(F, dat, dims)                        % dat is a simple Matlab ND array
        update_operation_(O, data, dims, old_data_op, varargin)    % data is an xplr.XData object
    end
    methods (Access='protected')
        function accepts_input(O, header)
            % Input header must match O.header_in for operation to apply.
            % This method can be overwritten in sub-classes for more
            % flexible acceptance of some differences.
            if ~isequal(header, O.header_in) % works also with non-scalar O
                error 'data header does not match operation specification'
            end
        end
        function b = change_dimension_id(O)
            % Whether output dimension header is intrinsically different
            % from input header, i.e. whether it corresponds to "something
            % else".
            % For example 2D ROI filtering is necessarily a dimension
            % change, but 1D ROI filtering isn't (both input and output
            % have the same label, lie in the same space, etc.). Performing
            % an FFT would be a dimension change even though the number of
            % dimensions is the same.
            % We consider that there is a dimension change when the number
            % of dimensions or the label(s) have changed.
            b = (O.nd_out ~= O.nd_in);
        end
    end
    methods
        function dim_id_out = get_dim_id_out(O, dim_id_in)
            % function dim_id_out = get_dim_id_out(O,dim_id_in)
            %---
            % generate an identifier for the replacing dimension(s) that
            % will be created when applying the operation to some
            % dimensions (identified by dim_id_in) of an xdata object.
            if O.change_dimension_id()
                dim_id_out = mod(sum(dim_id_in) + O.id_graph_node + (0:O.nd_out-1)*pi, 1);
            else
                dim_id_out = dim_id_in;
            end
        end
        function data = operation(O, data, dim_ids)
            % dimension number
            dims = data.dimension_number(dim_ids);
            % check input
            O.accepts_input(data.header(dims))            
            % actual code of operation will be in child class
            dat = data.data;                 % Matlab ND array
            dat = O.operation_(dat, dims);    % Matlab ND array
            % output header
            dim_id_out = O.get_dim_id_out(dim_ids);
            head = data.header;
            head(dims) = [];
            head = [head(1:min(dims) - 1), xplr.DimHeader(O.header_out, dim_id_out), head(min(dims):end)];
            % build output xdata object
            data = xplr.XData(dat, head);
        end
        function update_operation(O, data, dim_ids, old_data_op, varargin)
            % dimension number
            dims = data.dimension_number(dim_ids);
            % check input
            O.accepts_input(data.header(dims))            
            % actual code of operation will be in child class
            update_operation_(O, data, old_data_op, varargin{:});
        end
    end
    
    
    % Additional information in output header
    methods
        function [head_value, affected_columns] = set_add_header_info(F, head_value, add_header_info)
            affected_columns = [];
            for i=1:size(add_header_info, 2)
                label = add_header_info{1, i};
                values = add_header_info{2, i};
                if ~iscell(values)
                    if ischar(values), values = cellstr(values); else, values = num2cell(values); end
                end
                if ~isvector(values) || (~isscalar(values) && length(values) ~= size(head_value, 1))
                    error 'size of additional header info does not match number of selections'
                end
                
                % new label?
                idx = find(strcmp(label, {F.header_out.sub_labels.label}), 1);
                if isempty(idx)
                    % create new label
                    label_type = xplr.DimensionLabel.infer_type(values{1});
                    F.header_out = F.header_out.add_label(xplr.DimensionLabel(label, label_type));
                    idx = find(strcmp(label, {F.header_out.sub_labels.label}), 1);
                end
                
                % assign values
                head_value(:, idx) = values;
                affected_columns(end+1) = idx; %#ok<AGROW>
            end
        end
        function augment_header(F, new_label, label_type)
            if any(strcmp(new_label,{F.header_out.sub_labels.label})), return, end
            F.header_out = F.header_out.add_label(xplr.DimensionLabel(new_label, label_type));
        end
    end
    
    % Synchronization of operation definition in real world coordinates
    % system. (world_operation is the 'operation' property of a
    % worldOperand object)
    methods (Abstract)
        world_op = operation_data_to_space(O)       % get world operation based on opeartion definition in O
        update_operation_data_to_space(O, WO, event)   % updates WO.operation based on operation definition in O and argument event; must take care of launching WO 'ChangedOperation' event
        update_operation_space_to_data(O, world_operation, event)   % updates operation definition in O based on world operation and optional argument event
    end
    
    % Load/save
    methods (Abstract)
        copy_in(O, obj)   % copy the operation specification from another objec
    end
    methods
        function save_to_file(O, f_name)
            % function save_to_file(O,f_name)
            %---
            % save dataOperand object from file
            fn_savevar(f_name, O);
        end
        function load_from_file(O, f_name)
            % function load_from_file(O,f_name)
            %---
            % set current dataOperand object properties from information
            % saved in file (note that this does not replace object O, nor
            % affects any of the listener attached to it)
            
            % load from file
            obj = fn_loadvar(f_name);
            
            % checks
            if ~isa(obj, class(O))
                error('attempted to load a %s object, but file content is a %s', class(O), class(obj))
            end
            if ~isequal(obj.header_in, O.header_in)
                if ~isequal({obj.header_in.label}, {O.header_in.label})
                    error('operand loaded from file applies to dimensions %s, expected %s instead', fn_strcat({obj.header_in.label}, ','), fn_strcat({O.header_in.label}, ','))
                elseif ~isequal([obj.header_in.n], [O.header_in.n])
                    error('operand loaded from file applies on data of size %s, expected %s instead', num2str([obj.header_in.n], '%i '), num2str([O.header_in.n], '%i '))
                else
                    error('operand loaded from file does not apply to the same type of input headers as the current object')
                end
            end
            
            % copy property values
            O.copy_in(obj);
        end
    end
    
    % Context menu
    methods
        function context_menu(O, m)
            % function context_menu(O, m)
            %---
            % populate a context menu with actions that can be applied to
            % the filter
            % this function should be overwritten by sub-classes
            delete(get(m, 'children'))
            uimenu(m, 'enable', 'off', 'label', '(empty menu)')
        end
    end
    
end
