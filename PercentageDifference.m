function percentageDiff = PercentageDifference(A, B)
%PERCENTAGEDIFFERENCE
%   Function takes in two matrices and calculates the perecntage difference 
percentageDiff = (abs(A-B)./ ((A+B)./2)).*100;
end

