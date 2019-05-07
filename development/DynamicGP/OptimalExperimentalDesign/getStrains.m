function strains = getStrains(damage, load)
    P = py.sys.path;
    if count(P,'methods/2d_toy_problem/') == 0
        insert(P,int32(0),'methods/2d_toy_problem/');
    end
    % Input: 
    % .    N x 1 or 1 x N array of damage values between 0 and 1
    % .    N x 1 or 1 x N array of load values between 1 and 6 
    % Output: 
    % .    N x 1 array of strains
    N = length(damage);
    strains = zeros(N, 1);
    for idx = 1:N
        % note: 25 is hard coded, but you can change it.
        % it is the length of the beam.  
        strains(idx) = py.mymod2.pass_strains_to_matlab( ...
            damage(idx), load(idx), 25);
    end     
end