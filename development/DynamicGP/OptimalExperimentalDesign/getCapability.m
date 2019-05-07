function capability = getCapability(damage)
    P = py.sys.path;
    if count(P,'methods/2d_toy_problem/') == 0
        insert(P,int32(0),'methods/2d_toy_problem/');
    end
    % Input:
    % .    N x 1 or 1 X N array of damage values between 0 and 1
    N = length(damage);
    capability = zeros(N, 1);
    for idx = 1:N
        % note: 1, 25 are arbitrary, and useless.  I just
        % coded it poorly and I don't want to break my old
        % codes...
        capability(idx) = py.mymod2.pass_flutter_to_matlab(...
            damage(idx), 1, 25);
    end
end