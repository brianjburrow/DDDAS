function betas = invertMap(samps, refSamps, gammas, multi_indices)
    tarSamples = construct_vandermond(refSamps, ones(size(gammas)), multi_indices);
    betas = tarSamples \ samps(:,end);
end