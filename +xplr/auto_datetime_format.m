function format = auto_datetime_format(value_start, value_stop, min_step)
% function format = auto_datetime_format(value_start, value_stop, min_step)

% convert duration to date
is_duration = isduration(value_start);
if is_duration
    origin = datetime(0,0,0);
    value_start = origin + value_start;
    value_stop = origin + value_stop;
end

% date format: show day / year only if difference between start
% and end
if dateshift(value_start, 'start', 'day') == dateshift(value_stop, 'start', 'day')
    date_format = '';
elseif dateshift(value_start, 'start', 'year') == dateshift(value_stop, 'start', 'year')
    if is_duration
        date_format = 'd';
    else
        date_format = 'eee dd/MM';
    end
else
    if is_duration
        date_format = 'y';
    else
        date_format = 'eee dd/MM/yy';
    end
end
% if min_step is a duration, determine whether the "round"
% number should be calculated in seconds, minutes, hours or
% days
if min_step > years(.5)
    if isempty(date_format)
        if is_duration
            date_format = 'y';
        else
            date_format = 'dd/MM'; 
        end
    end
    time_format = '';
elseif min_step > days(.5)
    if isempty(date_format), date_format = 'dd/MM'; end
    time_format = '';
elseif min_step > hours(.5)
    time_format = 'hh:mm';
elseif min_step > minutes(.5)
    time_format = 'hh:mm';
else
    time_format = 'hh:mm:ss';
end
if isempty(date_format)
    format = time_format;
elseif isempty(time_format)
    format = date_format;
else
    format = [date_format ' ' time_format];
end

