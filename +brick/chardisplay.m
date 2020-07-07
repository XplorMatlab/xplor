function [out, errormsg] = chardisplay(varargin)
% function [str, errormsg] = chardisplay(val[,'ignorelimit'])
% function val = chardisplay(str,type)
%---
% build a string representation of a non-character value and converts back
% the string to the original type
% 
% Example: ...
% 
% See also brick.idx2str

% Thomas Deneux
% Copyright 2013-2017

if nargin==0, help brick.chardisplay, return, end

% Input
switch nargin
    case 1
        conversiondir = 'any->char';
        ignorelimit = false;
    case 2
        if strcmp(varargin{2},'ignorelimit')
            conversiondir = 'any->char';
            ignorelimit = true;
        else
            conversiondir = 'char->any';
        end
    otherwise
        error 'two many arguments'
end

% Conversion
errormsg = ''; 
switch conversiondir
    case 'any->char'
        % conversion anything -> character
        val = varargin{1};
        if ischar(val)
            str = val;
        elseif isnumeric(val)
            if ~islogical(val) && size(val,1)==1 && all(val>=0 & mod(val,1)==0) && all(diff(val)>0)
                % vector of increasing non-negative integers: try to arrange them smartly
                str = brick.idx2str(val);
            elseif (ndims(val)<=2 && numel(val)<=22) || ignorelimit
                str = num2str(val);
                str(:,end+1) = ';'; str(:,end+1) = ' ';
                str = reshape(str',1,numel(str));
                str(end-1:end) = [];
                str = regexprep(str,' *',' ');
            else
                str = [];
                errormsg = 'cannot display array with more than 2 dimensions or 22 elements';
            end
        else
            str = [];
            errormsg = sprintf('class ''%s'' cannot be represented has a string',class(val));
        end
        out = str;
    case 'char->any'
        % conversion character -> specified type
        str = varargin{1};
        type = varargin{2};
        switch type
            case 'char'
                val = str;
            case {'logical' 'double' 'single' 'int8' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64'}
                try
                    val = eval(['[' str ']']);
                    val = cast(val,type);
                catch %#ok<CTCH>
                    errormsg = 'string could not be evaluated';
                end
            otherwise
                errormsg = 'a string cannot be converted to class ''%s''';
        end
        out = val;
end

% Handle error
if ~isempty(errormsg) && nargout<2
    error(errormsg)
end
