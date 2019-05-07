classdef simple_benchmark
    
    properties
        none
    end
    
    methods(Static)
        function obj = simple_benchmark()
            obj.none = 0;
        end  
        function state = StateTransition(currentState, time)
                state = currentState + 2 + normrnd(0, 10^0.5);
        end
        function state = StateTransitionC(currentState, time)
                state = currentState + 2 ;
        end
        function data = GetData(currentState)
            data = 0.05 * currentState + ...
                normrnd(0, 2^0.5);
        end
        function data = GetDataC(currentState)
            data = 0.05 * currentState;
        end
    end
end