% Calculates the vector |\max_{i\ne x} \mu_i - \mu_x| over x.
function delta = Delta(mu)
	if length(mu) <= 1
		error('The argument to Delta() must be an array with more than one element');
    end
    M = length(mu);
    max_without = ones(size(mu));
    max_without(1) = max(mu(2:M));
    max_without(M) = max(mu(1:M-1));
	for x=[2:M-1]
		max_without(x) = max(max(mu(1:x-1)),max(mu(x+1:length(mu))));
    end
    delta = max_without - mu;
end
