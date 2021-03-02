% Error of 1*10e-4 relative to peak value 

function kielkopfApprox = Kielkopf(isoSize,dataSize,gammaL,gammaG,v,v0,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength,P,pShift)
[kiel,kielkopfApprox,tempLineStrength,eta,betaL,betaG,y,x,xh,Gx,Lx,K0y,t] = deal(cell(1,isoSize));
a1 = 0.254829592; a2 = -0.284496736;
a3 = 1.421413741; a4 = -1.453152027;
a5 = 1.061405429; p=0.3275911;
epsilon = 0.0990;
for n = 1:isoSize
    betaL{n} = gammaL{n};
    betaG{n} = gammaG{n};
    y{n} = (betaL{n}.*sqrt(log(2)))./betaG{n}; 
    xh{n} = ( y{n}./2 ) .*(1+ epsilon.*log(2) + sqrt( (1-epsilon*log(2))^2 + 4.*(log(2)./y{n}.^2) ) ) ;
    eta{n} = (y{n}.*xh{n})./(1+(y{n}.*xh{n}));
    t{n} = 1./(1+p.*y{n});
    K0y{n} = a1*t{n} + a2*t{n}.^2 + a3*t{n}.^3 + a4*t{n}.^4 + a5*t{n}.^5;
end
for n = 1:isoSize
    for k = 1:dataSize(n)
         x{n}(k,:) = (2*sqrt(log(2))./betaG{n}(k)).*(v{n}(k,:)-v0{n}(k)')-(P.*pShift{n}(k));
         Gx{n}(k,:) = exp(-log(2).*( (x{n}(k,:)) ./xh{n}(k)).^2);
         Lx{n}(k,:) = 1./(1+((x{n}(k,:))./xh{n}(k)).^2);
         kiel{n}(k,:) = (K0y{n}(k).*((1-eta{n}(k)).*Gx{n}(k,:) +eta{n}(k).*Lx{n}(k,:)));
         tempLineStrength{n}(k) = S_t0{n}(k) .*( (Q_tref(n)/Q_t(n)) .* (exp(-c2.*E_lower{n}(k)./T) ./ exp(-c2.*E_lower{n}(k)./T0))...
        .* ( (1-exp(-c2.*v0{n}(k)./T)) ./(1-exp(-c2.*v0{n}(k)./T0))));
         kielkopfApprox{n}(k,:) =  2*P*concentration*pLength./betaG{n}(k).*tempLineStrength{n}(k).*sqrt(log(2)/pi).*(kiel{n}(k,:));
    end
end

end

