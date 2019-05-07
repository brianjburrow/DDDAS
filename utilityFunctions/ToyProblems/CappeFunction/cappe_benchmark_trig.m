classdef cappe_benchmark_trig
    
    properties
        none
    end
    
    methods(Static)
        function obj = cappe_benchmark_trig()
            obj.none = 0;
        end  
        function state = StateTransition(currentState, time)
                state = currentState./2 + ...
                    15*currentState./(1 + currentState.^2) +...
                    4*cos(1.2*time) + normrnd(0, 0.5^0.5, 1, length(currentState));
        end
        function state = StateTransitionC(currentState, time)
                state = currentState./2 + ...
                    15*currentState./(1 + currentState.^2) +...
                    4*cos(1.2*time);
        end
        function data = GetData(currentState)
            data = 0.05 * currentState.^2 + 95*sin(1.5*currentState)/100 +...
                normrnd(0, 0.25, 1, length(currentState));
        end
        function data = GetDataC(currentState)
            data = 0.05 * currentState.^2 + 95*sin(1.5*currentState)/100;
        end
    end
end