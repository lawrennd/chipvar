function [data,originalC, originalB, originalGamma, originalAlpha,originalBeta]=chipVarFakeData(X,npts);
% CHIPVARFAKEDATA artificial data for checking chipVar

% CHIPVAR
randn('seed',73)
originalAlpha=0.85;
originalBeta=8.63;
nTrans=size(X,2);
originalGamma=5*randn(size(X,2),1);
factors=(cos(originalGamma)+ones(nTrans,1))/(2+4e-6)+1e-6*ones(nTrans,1);
originalB=randn(size(X))/sqrt(originalAlpha);
originalC=randn(nTrans,1);
for i=2:npts
    originalC=[originalC,sqrt(ones(size(X,2),1)-factors.^2).*...
        (randn(size(X,2),1))+factors.*originalC(:,end)];
end
data=(X.*originalB)*originalC+randn(size(X,1),npts)/sqrt(originalBeta);
       
        