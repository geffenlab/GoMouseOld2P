function [P] = binoCompare(xCorr,xTot,yCorr,yTot)

px = xCorr/xTot;
py = yCorr/yTot;

pcomb = (xCorr+yCorr)/(xTot+yTot);

Z =  (px-py)/sqrt(pcomb*(1-pcomb)*(1/xTot+1/yTot));

P = (1 - normcdf( abs(Z) ))/2; % two-tailed test on the Z value

% if P<0.05 % 
%     H=1;
% else
%     H=0;
% end

% Confidence interval 95%

SE = sqrt(((px*(1-px))/xTot)+((py*(1-py))/yTot));

ci(1) = px-py+(1.96*SE);
ci(2) = px-py-(1.96*SE);