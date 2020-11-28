
% Function models Gaussian probability density function
% given the set x, the mean of the data and the standard deviation 
function gaus = GaussianPDF(x,mean,stdDev)
gaus = (1/(stdDev*(sqrt(2*pi))))*exp((-1*(x-mean).^2/(2*stdDev.^2)));
end

