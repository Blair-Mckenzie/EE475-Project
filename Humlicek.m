function Kxy = Humlicek(x,y,x_k,sigma,a_k,b_k)
%HUMLICEK Summary of this function goes here
%   Detailed explanation goes here
Kxy = exp(-x.^2) + ( y./ ((x-x_k).^2 + sigma.^2) ) *( (b_k*( (x-x_k).^2 -sigma*(y + sigma)) - a_k*(x-x_k)*(y+2*sigma)) ./ ( (x-x_k).^2 + (y+sigma).^2)) ; 
end

