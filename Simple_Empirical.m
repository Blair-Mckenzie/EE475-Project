% Max width error of 0.55% at d = 0.1
% Error is 0.0118% or 0.01% at d = -+0.54

function simpleApprox = Simple_Empirical(isoSize,dataSize,gammaL,...
    gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P)
[approx,simpleApprox,tempLineStrength,d,c_L,c_G,gammaV,sigmaL,...
    sigmaG] = deal(cell(1,isoSize));
for n = 1:isoSize
  sigmaL{n} = gammaL{n}./2;
  sigmaG{n} = gammaG{n}./2;
  % Calculation of Voigt HMHW
  gammaV{n}=(0.5346.*gammaL{n})+sqrt(0.2166.*(gammaL{n}.^2)+gammaG{n}.^2);
  gammaV{n} = gammaV{n}./2;
  % Parameters required to calculate simple empirical approximation
  d{n} = ((sigmaL{n} - sigmaG{n})./(sigmaL{n} + sigmaG{n}));
  c_L{n}=0.68188+(0.61293.*d{n})-(0.18384.*(d{n}.^2))-(0.11568.*(d{n}.^3));
  c_G{n}=0.32460-(0.61825.*d{n})+(0.17681.*(d{n}.^2))+(0.12109.*(d{n}.^3));
  for k = 1:dataSize(n)
      approx{n}(k,:) = (c_L{n}(k).* 1/pi).*(gammaV{n}(k)./...
          ((v{n}(k,:)-v0{n}(k)).^2 + gammaV{n}(k).^2))+c_G{n}(k)...
          .*(sqrt(log(2))./sqrt(pi).*gammaV{n}(k)).*exp((-1*log(2)...
          .*(v{n}(k,:)-v0{n}(k)).^2 )./(gammaV{n}(k).^2)) ;
      % Calculation of temperature dependent line strength    
      tempLineStrength{n}(k) = S_t0{n}(k) .*( (Q_tref(n)/Q_t(n)).*...
          (exp(-c2.*E_lower{n}(k)./T) ./ exp(-c2.*E_lower{n}(k)./T0))...
      .* ((1-exp(-c2.*v0{n}(k)./T)) ./(1-exp(-c2.*v0{n}(k)./T0))));
      % Final aborsption for simple empirical approximation
      simpleApprox{n}(k,:) = (P*concentration*pLength./sigmaG{n}(k).*...
      tempLineStrength{n}(k).*sqrt(log(2)/pi).*(approx{n}(k,:)))./18;
  end
end
end