
% Function models Lorentzian probability density function
% given the set x,x0 and the half width half maxima 
function lorentz = LorentzianPDF(x,x0,HWHM)
lorentz = 1./((pi*HWHM)*(1+(((x-x0)/HWHM).^2)));
end

