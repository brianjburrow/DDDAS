% function i = argmax(list)
% Randomly choose from among the argmax set.
function i = argmax(list)
	best = max(list);
	i = ChooseRandomElement(find(list == best));
end

