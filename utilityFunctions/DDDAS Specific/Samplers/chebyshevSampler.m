function samples = chebyshevSampler(nSamples, lb, ub)
    samples = betarnd(0.5, 0.5, [nSamples, 1]);
    
    lbVec = lb*ones([nSamples, 1]);
    ubVec = ub*ones([nSamples, 1]);
    samples = (ubVec - lbVec).*(samples) + lbVec;
end