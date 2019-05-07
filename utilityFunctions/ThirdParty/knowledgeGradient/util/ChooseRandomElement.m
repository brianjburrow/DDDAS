% Choose an element from a list uniformly at random. 
function i = ChooseRandomElement(list)
	n = length(list);
	if (n == 1)
		i = list(1);
	else
		j = ceil(rand() * n);
		i = list(j);
	end
end

