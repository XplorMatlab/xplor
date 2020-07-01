classdef WorldOperand < xplr.GraphNode
% Similarly to dataOperand, a worldOperand object defines an operation to
% perform on data, but it is not associated with specific data headers;
% rather it will link data headers corresponding to the same measure space,
% for example time data headers with different temporal resolutions and
% offsets.
% When the bank creates a dataOperand object applying on measure headers,
% it immediately looks for (and creates if necessary) a worldOperand object
% to link this dataOperand object to.
% 
% See also xplr.dataOperand, xplr.Bank

properties (SetAccess='private')
    type        % type of operation: the class of the dataOperand objects it is linked to
    space_id     % identifier of the measure space the object operate on
end
properties
    operation   % operation definition: will be defined by the dataOperand objects it is linked to
end
   
events
    ChangedOperation
end
    

methods
    function wo = WorldOperand(do)
        % wo and do are respectively the constructed worldOperand object
        % and a first dataOperand object it links to

        % information about the type of operation and the space it operates
        % on
        wo.type = class(do);
        wo.space_id = do.header_in.get_measure_space_id();
        % connect dataOperand and worldOperand together
        wo.add_listener_exclusive_pair(do, ...
            'ChangedOperation', @(u,e)do.update_operation_space_to_data(wo.operation,e), ...
            'ChangedOperation', @(u,e)do.update_operation_data_to_space(wo,e));
        do.world_operand = wo;
        % obtain world operation by running dataOperand method
        wo.operation = do.operation_data_to_space();
    end
    function connect_data_operand(wo,do)
        % check
        if ~strcmp(class(do), wo.type) || ~isequal(do.header_in.get_measure_space_id(), wo.space_id)
            error 'cannot connect dataOperand to worldOperand as they are not compatible'
        end
        % connect
        wo.add_listener_exclusive_pair(do, ...
            'ChangedOperation', @(u,e)do.update_operation_space_to_data(wo.operation,e), ...
            'ChangedOperation', @(u,e)do.update_operation_data_to_space(wo,e));
        do.world_operand = wo;
        % set data operation
        do.update_operation_space_to_data(wo.operation)
    end
end

end
