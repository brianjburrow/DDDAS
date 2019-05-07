% Pass column vectors x and y.
function [variohat d sq] = EmpiricalVariogram(x, y)
assert(length(x) == length(y));
n = length(x);
d = zeros(1,n^2);
sq = zeros(1,n^2);
for k=1:n^2
	[i j] = ind2sub([n n], k);
	d(k) = norm(x(i)-x(j));
	sq(k) = (y(i)-y(j))^2;
end
uniqd = unique(d);
variohat = zeros(size(uniqd));
for k=1:length(uniqd)
	variohat(k) = mean(sq(find(d==uniqd(k))));
end
