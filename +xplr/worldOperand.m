classdef worldOperand < xplr.graphnode
% Similarly to dataOperand, a worldOperand object defines an operation to
% perform on data, but it is not associated with specific data headers;
% rather it will link data headers corresponding to the same measure space,
% for example time data headers with different temporal resolutions and
% offsets.
% When the bank creates a dataOperand object applying on measure headers,
% it immediately looks for (and creates if necessary) a worldOperand object
% to link this dataOperand object to.
% 
% See also xplr.dataOperand, xplr.bank

properties (SetAccess='private')
    type        % type of operation: the class of the dataOperand objects it is linked to
    spaceID     % identifier of the measure space the object operate on
end
properties
    operation   % operation definition: will be defined by the dataOperand objects it is linked to
end
   
events
    ChangedOperation
end
    

methods
    function WO = worldOperand(DO)
        % WO and DO are respectively the constructed worldOperand object
        % and a first dataOperand object it links to

        % information about the type of operation and the space it operates
        % on
        WO.type = class(DO);
        WO.spaceID = DO.headerin.getMeasureSpaceID();
        % connect dataOperand and worldOperand together
        WO.addListenerExclusivePair(DO, ...
            'ChangedOperation',@(u,e)DO.updateOperationSpace2Data(WO.operation,e), ...
            'ChangedOperation',@(u,e)DO.updateOperationData2Space(WO,e));
        % obtain world operation by running dataOperand method
        WO.operation = DO.operationData2Space();
    end
    function connectDataOperand(WO,DO)
        % check
        if ~strcmp(class(DO),WO.type) || ~isequal(DO.headerin.getMeasureSpaceID(), WO.spaceID)
            error 'cannot connect dataOperand to worldOperand as they are not compatible'
        end
        % connect
        WO.addListenerExclusivePair(DO, ...
            'ChangedOperation',@(u,e)DO.updateOperationSpace2Data(WO.operation,e), ...
            'ChangedOperation',@(u,e)DO.updateOperationData2Space(WO,e));
        % set data operation
        DO.updateOperationSpace2Data(WO.operation)
    end
end

end
