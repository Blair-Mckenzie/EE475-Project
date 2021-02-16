function simpleApprox = Simple_Empirical(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P)
[approx,simpleApprox,tempLineStrength,d,c_L,c_G,gammaV,sigmaL,sigmaG] = deal(cell(1,isoSize));

for n = 1:isoSize
    sigmaL{n} = gammaL{n}./2;
    sigmaG{n} = gammaG{n} ./2 ;
%     gammaV{n} = 0.5.*(C1.*sigmaL{n}+ sqrt(C2.*sigmaL{n}.^2 + 4*C3.*sigmaG{n}.^2));
    gammaV{n} =  (0.5346.*sigmaL{n} + sqrt(0.2166.*sigmaL{n}.^2 + sigmaG{n}.^2));
    d{n} = (sigmaL{n} - sigmaG{n})./(sigmaL{n} + sigmaG{n});
    c_L{n} = 0.68188 + 0.61293.*d{n} - 0.18384.*d{n}.^2 -0.11568.*d{n}.^3;
    c_G{n} = 0.32460 - 0.61825 .*d{n} + 0.17681.*d{n}.^2 + 0.12109.*d{n}.^3;
end

for n = 1:isoSize
    for k = 1:dataSize(n)
        approx{n}(k,:) = ( (c_L{n}(k) .* 1/pi) .* (gammaV{n}(k)./(v{n}(k,:)-v0{n}(k).^2) + gammaV{n}(k).^2) ) + c_G{n}(k) .* (sqrt(log(2))./ sqrt(pi)...
        .* gammaV{n}(k) ) .* exp( (-log(2).*(v{n}(k,:)-v0{n}(k)).^2 ) ./ (gammaV{n}(k).^2) ) ;
        tempLineStrength{n}(k) = S_t0{n}(k) .*( (Q_tref(n)/Q_t(n)) .* (exp(-c2.*E_lower{n}(k)./T) ./ exp(-c2.*E_lower{n}(k)./T0))...
        .* ( (1-exp(-c2.*v0{n}(k)./T)) ./(1-exp(-c2.*v0{n}(k)./T0))));
        simpleApprox{n}(k,:) =  2*P*concentration*pLength./sigmaG{n}(k).*tempLineStrength{n}(k).*sqrt(log(2)/pi).*(approx{n}(k,:));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Voigt Lineshape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Voigt HWHM approximation
% gammaV =  (0.5346.*gammaL + sqrt(0.2166.*gammaL.^2+gammaG.^2)./2);
% %Lorentzian and Gaussian HWHMs
% sigmaL = gammaL./2;
% sigmaG = gammaG ./2 ;
% %Variables used for Empirical Approximation
% d = (sigmaL - sigmaG)./(sigmaL + sigmaG);
% c_L = 0.68188 + 0.61293.*d - 0.18384.*d.^2 -0.11568.*d.^3;
% c_G = 0.32460 - 0.61825 .*d + 0.17681.*d.^2 + 0.12109.*d.^3;


% y for humlicek voigt
% y = (gammaL.*sqrt(log(2)))./sigmaG;
% n =12;
% sigma = 1.5;

% x_k = [ 0.314240376254359 0.947788391240164 1.597682635152605 2.279507080501060 3.020637025120890 3.889724897869782];
% w_k = [ 0.5701352362625 0.2604923102642 0.5160798561588 0.3905390584629 0.8573687043588 0.2658551684356];
% x_k = [-3.889724897869782 -3.020637025120890 -2.279507080501060 -1.597682635152605 -0.947788391240164 -0.314240376254359 0.314240376254359 0.947788391240164 1.597682635152605 2.279507080501060 3.020637025120890 3.889724897869782];
% w_k = [-0.2658551684356 -0.8573687043588 -0.3905390584629 -0.5160798561588 -0.2604923102642 -0.5701352362625 0.5701352362625 0.2604923102642 0.5160798561588 0.3905390584629 0.8573687043588 0.2658551684356];
 
% Vxy = zeros(step,4,dataSize);
% tempLineStrength = zeros(dataSize,1);
% result = zeros(step,n/2,dataSize);

    
    %     % x for humlicek voigt 
    %      x(k,:) = (2*sqrt(log(2))./sigmaG(k)).*(v(k,:)-v0(k)')-(P.*pShift(k));
    %      
    %     %empirical expression to approximate the Voigt function 
    %     simpleApprox(k,:) = ( (c_L(k) .* 1/pi) .* (gammaV(k)./(v(k,:)-v0(k).^2) + gammaV(k).^2) ) + c_G(k) .* (sqrt(log(2))./ sqrt(pi) .* gammaV(k) ) .* exp( (-log(2).*(v(k,:)-v0(k)).^2 ) ./ (gammaV(k).^2) ) ;
    %     
    %     vf = real(fadf(x(k,:)+ 1i.*y(k)));
    %     
    %     for kk = 1:n
    %         a_k = -1*(1/pi)*w_k(kk)*exp(sigma^2)*sin(2*x_k(kk)*sigma);
    %         b_k = (1/pi)*w_k(kk)*exp(sigma^2) * cos(2*x_k(kk)*sigma);
    %         result(:,kk,k) =  ( y(k)./ ((x(k,:)-x_k(kk)).^2 + sigma^2) ) .*( (b_k.*( ((x(k,:)-x_k(kk)).^2) -(sigma.*(y(k) + sigma))) - (a_k.*(x(k,:)-x_k(kk)).*(y(k)+2.*sigma)) ) ./ ( (x(k,:)-x_k(kk)).^2 + (y(k)+sigma).^2)) ;
    %     end
    %     humlicekApprox(k,:) = sum(result(:,:,k)');
    %     humlicekApprox(k,:) = exp(-1.*x(k,:).^2) +  sum(result(:,:,k)');


end