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
    
    % Constructor
    methods
        
    end
    
    % Size
    methods
        function x = get.reductionfactor(O)
            x = prod([O.headerin.n])/prod([O.headerout.n]);
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
    end
    
    % Operation
    methods (Abstract, Access='protected')
        dat = operation_(F,dat,dims)                        % dat is a simple Matlab ND array
        updateOperation_(O,data,dims,olddataop,varargin)    % data is an xplr.xdata object
    end
    methods (Access='protected')
        function checkdata(O,data,dims)
            % input header must match O.headerin! 
            if ~isequal(data.header(dims),O.headerin) % works also with non-scalar O
                error 'data header does not match operation specification'
            end
        end
        function b = changedimensionID(O)
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
            b = (O.ndout ~= O.ndin) || ~all(strcmp({O.headerout.label}, {O.headerin.label}));
        end
    end
    methods
        function dimIDout = getdimIDout(O,dimIDin)
            % function dimIDout = getdimIDout(O,dimIDin)
            %---
            % generate an identifier for the replacing dimension(s) that
            % will be created when applying the operation to some
            % dimensions (identified by dimIDin) of an xdata object.
            if O.changedimensionID()
                dimIDout = mod(sum(dimIDin) + O.idGraphNode + (0:O.ndout-1)*pi, 1);
            else
                dimIDout = dimIDin;
            end
        end
        function data = operation(O,data,dimIDs)
            % dimension number
            dims = data.dimensionNumber(dimIDs);
            % check input
            checkdata(O,data,dims)            
            % actual code of operation will be in child class
            dat = data.data;                 % Matlab ND array
            dat = O.operation_(dat,dims);    % Matlab ND array
            % output header
            dimIDout = O.getdimIDout(dimIDs);
            head = data.header;
            head(dims) = [];
            head = [head(1:dims(1)-1) xplr.dimheader(O.headerout,dimIDout) head(dims(1):end)];
            % build output xdata object
            data = xplr.xdata(dat,head);
        end
        function updateOperation(O,data,dimIDs,olddataop,varargin)
            % dimension number
            dims = x.dimensionNumber(dimIDs);
            % check input
            checkdata(F,x,dims)            
            % actual code of operation will be in child class
            updateOperation_(O,data,olddataop,varargin{:});
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