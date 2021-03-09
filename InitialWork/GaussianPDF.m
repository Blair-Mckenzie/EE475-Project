
% Function models Gaussian probability density function
% given the set x, the mean of the data and the standard deviation 
function gaus = GaussianPDF(x,mean,alpha)
gaus = (1/(alpha*(sqrt(2*pi))))*exp((-1*(x-mean).^2/(2*alpha.^2)));
% gaus = (exp(-1.*(x.^2)./2*alpha^2))./ (alpha *sqrt(2*pi)) ;
% gaus = sqrt(log(2) / pi) / alpha * exp(-(x ./ alpha).^2 .*log(2));
end

