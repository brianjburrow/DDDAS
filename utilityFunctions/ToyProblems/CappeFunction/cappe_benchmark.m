classdef cappe_benchmark
    
    properties
        none
    end
    
    methods(Static)
        function obj = cappe_benchmark()
            obj.none = 0;
        end  
        function state = StateTransition(currentState, time)
                state = currentState./2 + ...
                    25*currentState./(1 + currentState.^2) +...
                    8*cos(1.2*time) + normrnd(0, 10^0.5, 1, length(currentState));
        end
        function state = StateTransitionC(currentState, time)
                state = currentState./2 + ...
                    25*currentState./(1 + currentState.^2) +...
                    8*cos(1.2*time);
        end
        function data = GetData(currentState)
            data = 0.05 * currentState.^2 + ...
                normrnd(0, 1, 1, length(currentState));
        end
        function data = GetDataC(currentState)
            data = 0.05 * currentState.^2;
        end
    end
end