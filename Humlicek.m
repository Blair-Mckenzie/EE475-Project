function Kxy = Humlicek(x_k,sigma,w_k)

a_k = -1*(1/pi)*w_k*exp(sigma)^2 * sin(2*x_k)*sigma;
b_k = (1/pi)*w_k*exp(sigma)^2 * cos(2*x_k)*sigma ;
n =12;
sum =0;
Kxy = exp(-x.^2) + ( y./ ((x-x_k).^2 + sigma.^2) ) *( (b_k*( (x-x_k).^2 -sigma*(y + sigma)) - a_k*(x-x_k)*(y+2*sigma)) ./ ( (x-x_k).^2 + (y+sigma).^2)) ; 
for k = -n/2:1:n/2

    
end


end

