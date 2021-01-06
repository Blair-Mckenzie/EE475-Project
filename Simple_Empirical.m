function simpleApprox = Simple_Empirical(t,Y,X)


K_function = ((exp(-t.^2))./ ((Y.^2) + (X - t).^2));
Kxy = Y/pi * int(K_function,t,-Inf,Inf);

L_function = ((X-t) .* exp(-t.^2))./ (Y.^2 + (X-t).^2);
Lxy = 1/pi * integral(L_function,-Inf,Inf);



end