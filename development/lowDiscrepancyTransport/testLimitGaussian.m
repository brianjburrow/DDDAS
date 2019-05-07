clear all 
close all
clc

ss = sobolset(2);
samples = ss(1:20000000,:);

samples = uniformTOgaussian(samples);

scatter(samples(:,1), samples(:,2), 10, 'black', 'filled')

function stdNormSample = uniformTOgaussian(uniSample)
    stdNormSample = norminv(uniSample);
end