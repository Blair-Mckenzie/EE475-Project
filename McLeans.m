% Works reasonably well for large Y
% For y = 1 relative error is less than 1*10e-2
% For y = 0.001 the relative eror reaches almost 10

function mcleans = Mcleans(isoSize,dataSize,gammaG,v,v0,P,pShift,Y,S_t0,Q_tref,Q_t,c2,E_lower,T,T0,concentration,pLength)
%Mclean's Coefficients
A = [-1.215, -1.3509, -1.215, -1.3509];
B = [1.2359, 0.3786, -1.2359, -0.3786];
C = [-0.3085, 0.5906, -0.3085, 0.5906];
D = [0.021, -1.1858, -0.021, 1.1858];

[X,mcleans,Vxy,tempLineStrength] = deal(cell(1,isoSize));

for n = 1:isoSize
    for k = 1:dataSize(n)
        %Calculating X for Voigt lineshape
        X{n}(k,:) = (2*sqrt(log(2))./gammaG{n}(k)).*(v{n}(k,:)-v0{n}(k)')-(P.*pShift{n}(k));  
        
        for index = 1:4
            Vxy{n}(:,index,k) = ((C(index).*(Y{n}(k)-A(index)))+D(index).*(X{n}(k,:)-B(index))) ./ ((Y{n}(k)-A(index)).^2 + (X{n}(k,:)-B(index)).^2);
        end
        tempLineStrength{n}(k) = S_t0{n}(k) .*( (Q_tref(n)/Q_t(n)) .* (exp(-c2.*E_lower{n}(k)./T) ./ exp(-c2.*E_lower{n}(k)./T0))...
        .* ( (1-exp(-c2.*v0{n}(k)./T)) ./(1-exp(-c2.*v0{n}(k)./T0))));
        mcleans{n}(k,:) =  2*P*concentration*pLength./gammaG{n}(k).*tempLineStrength{n}(k).*sqrt(log(2)/pi).*sum(Vxy{n}(:,:,k).');
    end
end
    
end