%% Template for using transport_map_mcmc class
% Author: Brian J. Burrows, Texas A&M University 2018

rng(0)                                                                      % Fix random number generator if desired
measurement = 2;                                                            % Fixed measurement (mcmc_samples.mat depends on this measurement)
%% Define MCMC properties
mcmcProp.adaptGap = 10;                                                     % Number of MCMC samples generated between Transport map Updates
mcmcProp.numTotalSamples = 20000;                                           % Scalar value, total number of MCMC samples desired
imported = load("mcmc_samples.mat");
imported.samples = imported.samples(1:5000,:);

sortIndices = determine_map_ordering(imported.samples);                     % MUST REORDER LOG_PDF FUNCTIONS, FORWARD MODEL, ETC. IF SORTING SAMPLES.  May be necessary if dimension > 1

mcmcProp.samples = imported.samples + 50*ones(size(imported.samples)) ;     % [numSamples x numDimensions], numSamples << numTotalSamples

%% Define Transport Map Properties
tmProp.multi_index_type = "TO";                                             % Multi index set type: "TO", "NM", "NC" currently supported
tmProp.polynomialOrder = 4;                                                 % Maximum hermite polynomial order used in expansion.  Max 4 supported.  5+ Gives weird results

opt = optimoptions(...                                                      % Set optimizer options.
                'fmincon', 'MaxIterations',  1000, 'algorithm',...
                'sqp', 'StepTolerance', 1000*eps,...
                'MaxFunctionEvaluations', 10000,...
                'SpecifyObjectiveGradient', true,...
                'SpecifyConstraintGradient',true,...
                'UseParallel',false);                                       % Parallel not useful since we use analytic gradient information...
            
%% Define Target Pdf properties
tarPdf.log_target_pdf = @(xx) log_likelihood(measurement, xx) + log_prior(xx);      


%%%%%%%%%%%%%%%%%%% End User Input %%%%%%%%%%%%%%%%%%%%
%% Run Transport Map Accelerated MCMC
object = transport_map_mcmc(mcmcProp, tmProp, opt, tarPdf);

%% Plot Sampling Results
object.plotMarginalHists()

object.plotTrace()

object.plotReferenceDistributions()                                         % First subplot is target, Second should look like gaussian, third should look like target distribution
function p = log_likelihood(data, state)
    mu = forwardModel(state);
    p  = logmvnpdf(data, mu, 1); %+ logmvnpdf(state, 2, 1);
end

function log_p = log_prior(state)
    log_p = logmvnpdf(state, 50, 10);
end

function meas = forwardModel(input)
    %meas = 1*input + sin(input);
    meas = (input-50)^2;
end