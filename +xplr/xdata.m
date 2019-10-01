classdef xdata < xplr.graphnode
    % function x = xdata(dat[,head[,name]])
    %---
    % A container for data associated with header information
    % A number of different events are thrown when the data or header
    % information is being changed    
    % When headers are not provided, opens a graphic interface allowing
    % user to set the headers.
    % 
    % Input:
    % dat   ND array
    % head  a cell array with as many elements as data dimensions, 
    %       each element is itself a cell array containing arguments for
    %       the xplr.header constructor
    % name  string
    %
    % See also xplr.header
    
    properties (SetAccess='private')
        data        % ND array
        header = xplr.header.empty(1,0);
        name = '';
    end
    properties (Dependent, SetAccess='private')
        sz
        nd
    end
    
    events
        ChangedData % sent with info xplr.eventinfo('data',chghead)
    end
    
    methods
        function x = xdata(dat,head,name)
            if nargin==0, return, end % default empty data            
            if nargin<2 || isempty(head)
                % open user edition window to edit headers
                head = xplr.editHeader(dat);
            elseif ischar(head)
                head = {head};
            end
            if iscell(head)
                % create xplr.header objects from labels
                labels = head;
                head = xplr.header.empty(1,0);
                for i=1:length(labels)
                    if iscell(labels{i})
                        head(i) = xplr.header(labels{i}{:});
                    else
                        head(i) = xplr.header(labels{i},size(dat,i));
                    end
                end
            end
            x.updateDataDim('global',[],dat,head)
            if nargin>=3
                x.name = name;
            end
        end
        function y = copy(x)
            y = xplr.xdata(x.data,x.header);
        end
        function nd = get.nd(x)
            nd = length(x.header);
        end
        function s = get.sz(x)
            s = [x.header.n];
        end
        function chgData(x,data)
            if isequal(data,x.data), return, end
            % changes in size are allowed only in the 'measure' dimensions
            datasz = xplr.strictsize(data,x.nd);
            if length(datasz)>x.nd, error 'Cannot increase number of dimensions. Use updateData to change both data and headers.', end
            chgsz = (datasz~=x.sz);
            if any([x.header(chgsz).iscategoricalwithvalues]), error 'Cannot change data size in categorical dimensions. Use updateData to change both data and headers', end
            % set data
            if ~isreal(data), error 'data cannot be complex', end
            x.data = data;
            if ~any(chgsz)
                notify(x,'ChangedData',xplr.eventinfo('data','chgdata'))
                return
            end
            % update dimension headers if necessary
            for i=find(chgsz)
                x.header(i) = updateMeasureHeader(x.header(i),datasz(i));
            end
            notify(x,'ChangedData',xplr.eventinfo('data','chgdim',find(chgsz)))
        end
        function updateData(x,flag,dim,ind,value,newhead)
            % function updateData(x,flag,dim,ind,value,newhead)
            %---
            % arguments value and newhead are supposed to be only the
            % updated parts, i.e.:
            % - for flags 'new' and 'chg', value is the data only for
            %   indices ind
            % - newhead is only the header in dimension dim
            % however giving the full updated data for value, or the full
            % header for newhead, is tolerated
            
            % check that value is real
            if nargin>=5 && ~isreal(value), error 'data cannot be complex', end
            
            % update header
            if nargin>=6
                if length(newhead)==x.nd
                    % giving the full headers instead of only the updated one
                    newhead = newhead(dim);
                end
                checkHeaderUpdate(x.header(dim),flag,ind,newhead)
            else
                newhead = updateHeader(x.header(dim),flag,ind);
            end
            if strcmp(flag,'all') && isequal(newhead,x.header(dim))
                % flag 'chgdata' might be preferable to 'all' to indicate
                % that header did not change
                flag = 'chgdata'; 
            end
            tmp = newhead; newhead = x.header; newhead(dim) = tmp; clear tmp
            % update data
            if ~fn_ismemberstr(flag,{'chgdata' 'all' 'chgdim'})
                s = substruct('()',repmat({':'},1,x.nd));
                s.subs{dim} = ind;
            end
            switch flag
                case {'all' 'chgdata' 'chgdim'}
                    newdata = value;
                case {'chg' 'new' 'chg&new' 'chg&rm'}
                    if size(value,dim)==newhead(dim).n
                        % giving the full data instead of only the updated part
                        newdata = value;
                    else
                        if strcmp(flag,'chg&new')
                            s.subs{dim} = [ind{:}];
                        elseif strcmp(flag,'chg&rm')
                            s.subs{dim} = ind{1};
                        end
                        newdata = subsasgn(x.data,s,value);
                        if strcmp(flag,'chg&rm')
                            s.subs{dim} = ind{2};
                            newdata = subsasgn(newdata,s,[]);
                        end
                    end
                case 'remove'
                    newdata = subsasgn(x.data,s,[]);
                case 'perm'
                    newdata = subsref(x.data,s);
                otherwise
                    error('invalid flag ''%s'' for xdata updateData method')
            end
            if ~isequal(xplr.strictsize(newdata,length(newhead)),[newhead.n])
                error 'new data size does not match header(s)'
            end
            % really update data only now (after all checks occured)
            x.header = newhead;
            x.data = newdata;
            % notification
            notify(x,'ChangedData',xplr.eventinfo('data',flag,dim,ind))
        end
        function updateDataDim(x,flag,dim,newdata,newhead)
            % update header
            switch flag
                case 'global'
                    x.header = newhead;
                case 'chgdim'
                    if length(dim)~=length(newhead), error 'length of new header does not match number of new dimensions', end
                    x.header(dim) = newhead;
                case 'insertdim'
                    if length(dim)~=length(newhead), error 'length of new header does not match number of new dimensions', end
                    ndnew = x.nd + length(dim);
                    newhead = [x.header newhead];
                    perm = zeros(1,ndnew); 
                    perm(setdiff(1:x.nd,dim)) = 1:x.nd;
                    perm(dim) = x.nd + (1:length(dim));
                    x.nd = ndnew;
                    x.header = newhead(perm);
                case 'rmdim'
                    x.header(dim) = [];
                case 'permdim'
                    perm = dim;
                    if ~isequal(unique(perm),1:x.nd), error 'permutation does not match number of dimensions', end
                    x.header = x.header(perm);
                    if nargin>=5 && ~isequal(newhead,x.header), error 'new header does not match dim permutation', end
                otherwise
                    error('invalid flag ''%s'' for xdata updateDataDim method')
            end
            % update data
            if ~isreal(newdata), error 'data cannot be complex', end
            switch flag
                case 'permdim'
                    x.data = permute(x.data,perm);
                    if nargin>=4 && ~isequal(newdata,x.data), error 'new data does not match dim permutation', end
                otherwise
                    if x.nd==1 && isvector(newdata), newdata = newdata(:); end % don't generate an error for a row vector input
                    if ~isequal(xplr.strictsize(newdata,x.nd),x.sz)
                        error 'new data size does not match new header'
                    end
                    x.data = newdata;
            end
            % notification
            notify(x,'ChangedData',xplr.eventinfo('data',flag,dim))
        end
    end
    
end


