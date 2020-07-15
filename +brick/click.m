function click(hobj,value)
%CLICK Mimick the effect of clicking an object or hitting a key inside a given figure
%---
% function click(pushbutton|togglebutton)
% function click(edit|radiobutton|listbox|etc.,value)
% function click(figure|axis|etc.,newCurrentPoint|keyPressed)
%---
% Mimick the effect of user clicking (or hitting a key) on a specified
% object and raising the appropriate callback.

switch get(hobj,'type')

    case 'uicontrol'
        switch get(hobj,'style')
            case 'pushbutton'
                set(hobj,'value',1)
            case 'togglebutton'
                set(hobj,'value',~get(hobj,'value'))
            case 'edit'
                if ischar(value)
                    str = value;
                elseif isnumeric(value) || islogical(value)
                    str = num2str(value);
                else
                    error 'value must be char, numeric or logical'
                end
                set(hobj,'string',str)
            otherwise
                error('object ''%s'' not handled yet',get(hobj,'style'))
        end
        callback = get(hobj,'callback');
        brick.evalcallback(callback,hobj,[])
        
    case {'figure' 'uipanel' 'axis'}
        if isnumeric(value)
            if length(value)~=2
                error 'new CurrentPoint value must have 2 elements'
            end
            set(hobj,'CurrentPoint',value)
            callback = get(hobj,'buttondownfcn');
            brick.evalcallback(callback,hobj,[])
            if strcmp(type,'figure')
                callback = get(hobj,'WindowButtonDownFcn');
                brick.evalcallback(callback,hobj,[])
            end
        elseif ischar(value)
            set(hobj,'CurrentCharacter',value)
            e = struct('Character',value,'Modifier',{{}},'Key',lower(value), ...
                'Source',hobj,'EventName','KeyPress');
            if isscalar(value) && value~=lower(value)
                e.Modifier = {'shift'};
            end
            callback = get(hobj,'keypressfcn');
            brick.evalcallback(callback,hobj,e)
        end
        
    otherwise
        error('object ''%s'' not handled yet',get(hobj,'type'))
        
end
