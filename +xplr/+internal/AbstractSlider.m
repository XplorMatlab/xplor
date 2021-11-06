classdef AbstractSlider < handle
    % Functionality of the slider, but does not take care of its display
    % NOT USED - NOT CLEAR WHAT USAGE
    
    properties (Access=private)
        % internal representation with integer of a descretized log scale: 
        % -10 ->   0.1
        %   0 ->   1
        %   1 ->   1.25
        %   2 ->   1.6
        %  10 ->  10
        %  20 -> 100
        idx_min
        idx_max
        n_step
        idx
    end
    properties (Dependent)
        value_min
        value_max
        value       % value
        position    % graphic position: between 0 and 1
    end
    
    % Constructor
    methods
        function S = AbstractSlider(value_min, value_max, value)
            S.set_min_max(S.index2value(S.value_min), S.index2value(S.value_max))
            S.value = value;
        end
    end
    
    % Conversion between internal indices and value
    methods
        function idx = value2index(S, val)
            assert(val > 0)
            x = log10(val);
            idx = round(10 * x);
        end
        function val = index2value(S, idx)
            steps10 = [1, 1.25, 1.6, 2, 2.5, 3, 4, 5, 6, 8];
            val = 10.^floor(idx/10) .* steps10(mod(idx, 10)); 
        end
        function value = get.value_min(S)
            value = S.index2value(S.value_min);
        end
        function value = get.value_max(S)
            value = S.index2value(S.value_max);
        end
        function set_min_max(S, min, max)
            assert(min < max)
            S.idx_min = min;
            S.idx_max = max;
            S.n_step = max - min;
        end
        function set.value_min(S, val)
            S.set_min_max(S.value2index(val), S.idx_max)
        end
        function set.value_max(S, val)
            S.set_min_max(S.idx_min, S.value2index(val))
        end
        function val = get.value(S)
            val = S.index2value(S.idx);
        end
        function set.value(S, val)
            S.idx = fn_coerce(S.value2index(val), S.value_min, S.value_max);
        end
    end
    
    % Slider position
    methods
        function pos = get.position(S)
            pos = (S.idx - S.idx_min) / (S.idx_max - S.idx_min);
        end
        function set.position(S, pos)
            pos = fn_coerce(pos, 0, 1);
            S.idx = round(S.idx_min + pos * (S.idx_max - S.idx_min));
        end
        function value = step(S, n)
            S.idx = fn_coerce(S.idx + round(n), S.idx_min, S.idx_max);
            if nargout > 0
                value = S.value;
            end
        end
    end
end