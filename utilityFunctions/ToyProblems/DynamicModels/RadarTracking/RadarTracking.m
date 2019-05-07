classdef RadarTracking < handle
    %RADARTRACKING Summary of this class goes here
    %   Implementation of the radar tracking problem from 
    %   "karlgaard2006comparison.pdf included in this folder
    
    properties
        a                                                                   % km
        b                                                                   % km
        eta                                                                 % m^-1
        sigma_1                                                             % m, standard deviation of 1st gaussian mixture component
        sigma_2                                                             % m, standard deviation of 2nd gaussian mixture component
        gamma                                                               % unitless
        c_2                                                                 % unitless
        epsilon                                                             % unitless
        x                                                                   % [ km; km; 1/sqrt(m)]
        P                                                                   % State error covariace matrix
        maxTime                                                             % seconds
    end
    
    methods
        function obj = RadarTracking()
            %RADARTRACKING Construct an instance of this class
            %   Folowing initialization from karlgaard2006comparison.pdf
            obj.a = 30.5;
            obj.b = 30.5;
            obj.eta = 1.64*10^-4;
            obj.sigma_1 = 30.5;
            obj.sigma_2 = 5 * obj.sigma_1;
            obj.gamma = 1.345;
            obj.c_2 = 3.0;
            obj.x = [91.5;6.1;0.01];         
            obj.P = diag([0.31^2;0.06^2;0.02^2]);
            obj.maxTime = 60;          
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function outputArg = takeMeasure(obj, inputArg)
            % Take a noisy radar measurement of the system
            bin = binornd(1, 0.5);
            if bin == 0
                % sample from first mixture component
            else
                % sample from second mixture component
            end
            outputArg = sqrt(obj.b^2 + (inputArg(1) - obj.a).^2) + w_k;
        end
    end
end

