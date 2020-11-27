
% Function returna FWHM Gaussian
% Given the speed of light(c), Boltzmann constant(Kb), central frequency(v0),
% Temperature of the system(T) and the molecular mass(m)
function gammaG = GammaDoppler(v0,T,M)
gammaG = v0.*7.1623e-7*sqrt(T/M);
end
