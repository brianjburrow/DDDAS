% logy = LogNormPDF(z)
% Returns the log of the normal pdf evaluated at z.  z can be a vector or a scalar.
function logy = LogNormPDF(z)
	const = -.5*log(2*pi); % log of 1/sqrt(2pi).
	logy = const - z.^2/2;
end
