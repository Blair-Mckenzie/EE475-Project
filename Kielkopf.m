function kielkopfApprox = Kielkopf(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P,pShift)
[kiel,kielkopfApprox,tempLineStrength,eta,betaL,betaG,y,x,xh,Gx,Ex,Lx,vf] = deal(cell(1,isoSize));
% b0 = 0.47047; % b1 = 0.61686;% b2 =-0.16994;% b3 = 1.32554
% C1 = 1.0692;% C2 = 0.86639;% C3 = 1;
epsilon = 0.0990;
for n = 1:isoSize
    betaL{n} = gammaL{n}./2;
    betaG{n} = gammaG{n} ./2 ;
    y{n} = (betaL{n}.*sqrt(log(2)))./betaG{n}; 
    xh{n} = ( y{n}./2 ) .*(1+ epsilon.*log(2) + sqrt( (1-epsilon*log(2))^2 + 4.*(log(2)./y{n}.^2) ) ) ;
    eta{n} = (y{n}.*xh{n})./(1+(y{n}.*xh{n}));
    
end

for n = 1:isoSize
    for k = 1:dataSize(n)
         % x for humlicek voigt 
         x{n}(k,:) = (2*sqrt(log(2))./betaG{n}(k)).*(v{n}(k,:)-v0{n}(k)')-(P.*pShift{n}(k));
         vf{n} = real(fadf(x{n}(k,:)+ 1i.*y{n}(k)));
         Gx{n}(k,:) = exp(-log(2).*( (x{n}(k,:)) ./xh{n}(k)).^2);
         Lx{n}(k,:) = 1./(1+((x{n}(k,:))./xh{n}(k)).^2);
         Ex{n}(k,:) = (0.8029- 0.4207.*x{n}(k,:).^2) ./ (1+0.2030.*x{n}(k,:).^2+0.07335.*x{n}(k,:).^4) ; 
         kiel{n}(k,:) = (1-eta{n}(k)).*Gx{n}(k,:) +eta{n}(k).*Lx{n}(k,:)+eta{n}(k).*(1-eta{n}(k)).*Ex{n}(k,:).*(Gx{n}(k,:)-Lx{n}(k,:));
         tempLineStrength{n}(k) = S_t0{n}(k) .*( (Q_tref(n)/Q_t(n)) .* (exp(-c2.*E_lower{n}(k)./T) ./ exp(-c2.*E_lower{n}(k)./T0))...
        .* ( (1-exp(-c2.*v0{n}(k)./T)) ./(1-exp(-c2.*v0{n}(k)./T0))));
         kielkopfApprox{n}(k,:) =  P*concentration*pLength./betaG{n}(k).*tempLineStrength{n}(k).*sqrt(log(2)/pi).*(kiel{n}(k,:));
%         kielkopfApprox{n} =  real(fadf(x{n}(k,:)+ 1i.*y{n}(k)));
    end
end


end

