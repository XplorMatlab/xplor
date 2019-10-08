classdef dataoperand < xplr.graphnode
% dataoperand
% Abstract class defining an operation on an xdata object
    
    properties (SetAccess='protected')
        headerin
        headerout
    end
    properties (Dependent, SetAccess='private')
        szin
        ndin
        szout
        ndout
    end
    properties (Dependent, SetAccess='private')
        reductionfactor
    end
    % properties below are not handled by dataoperand class and
    % sub-classes, but rather by the objects that use them
    properties 
        linkkey = 0
        %shared = struct;
    end
   
    events
        ChangedOperation
    end
    
    methods (Abstract)
        dataop = operation(O,data,dims) % code should start with 'O.checkdata(data,dims)'
        updateOperation(O,data,dims,olddataop,evnt)
    end
    
    % Get/Set Dependent
    methods
        function x = get.reductionfactor(O)
            x = prod([O.headerin.n])/prod([O.headerout.n]);
        end
    end
    
    methods
        function checkdata(O,data,dims)
            % input header must match O.headerin! 
            if ~isequal(data.header(dims),O.headerin) % works also with non-scalar O
                error 'data header does not match operation specification'
            end
        end
        function sz = get.szin(O)
            sz = [O.headerin.n];
        end
        function nd = get.ndin(O)
            nd = length(O.headerin);
        end
        function sz = get.szout(O)
            sz = [O.headerout.n];
        end
        function nd = get.ndout(O)
            nd = length(O.headerout);
        end
        function dimbef2aft = followdims(O,ndbef,dims)
        % function dimbef2aft = followdims(O,nd,dims)
        % for each dimension of an original nd-dimensional data, at
        % which new position is it going in the operated data
        
            if length(dims)~=O.ndin, error 'number of dimensions does not match filter input header', end
            [otherdim1 otherdim2] = deal(1:dims(1)-1,setdiff(dims(1):ndbef,dims));
            dimbef2aft = zeros(1,ndbef);
            dimbef2aft(otherdim1) = otherdim1;
            switch O.ndout
                case 0
                    dimbef2aft(dims) = 0;
                case 1
                    dimbef2aft(dims) = dims(1);
                otherwise
                    error 'case not handled'
            end
            ndaft = ndbef + (O.ndout-O.ndin);
            dimbef2aft(otherdim2) = dims(1)+O.ndout:ndaft;
        end
    end
    
    % Additional information in output header
    methods
        function [headvalue affectedcolumns] = setAddHeaderInfo(F,headvalue,addheaderinfo)
            affectedcolumns = [];
            for i=1:size(addheaderinfo,2)
                label = addheaderinfo{1,i};
                values = addheaderinfo{2,i};
                if ~iscell(values), values = num2cell(values); end
                if ~isvector(values) || (~isscalar(values) && length(values)~=size(headvalue,1))
                    error 'size of additional header info does not match number of selections'
                end
                
                % new label?
                idx = find(strcmp(label,{F.headerout.sublabels.label}),1);
                if isempty(idx)
                    % create new label
                    labeltype = xplr.dimensionlabel.infertype(values{1});
                    F.headerout = addLabel(F.headerout,xplr.dimensionlabel(label,labeltype));
                    idx = find(strcmp(label,{F.headerout.sublabels.label}),1);
                end
                
                % assign values
                headvalue(:,idx) = values;
                affectedcolumns(end+1) = idx; %#ok<AGROW>
            end
        end
        function augmentHeader(F,newlabel,labeltype)
            if any(strcmp(newlabel,{F.headerout.sublabels.label})), return, end
            F.headerout = addLabel(F.headerout,xplr.dimensionlabel(newlabel,labeltype));
        end
    end
    
end