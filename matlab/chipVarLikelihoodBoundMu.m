function f=chipVarLikelihoodBoundMu(model,data,X,...
                      expectationsC, expectationsB,expectationsMu)
% CHIPVARLIKELIHOODBOUNDMU variational lower bound on marginal likelihood

% CHIPVAR
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
beta=model.beta;
alpha=model.alpha;
Gamma=model.Gamma;
factors=(cos(Gamma)+ones(nTrans,1))/(2+4e-6)+1e-6*ones(nTrans,1);
factors=0.99*factors;
%factors=(1-1e-6)*exp(Gamma)./(ones(nTrans,1)+exp(Gamma));
% Gamma(find(factors==1))=Gamma(find(factors==1))-1e-3;
% factors=cos(Gamma).^2;
%expectationsB.bChi=expectationsB.b.*X;
condBit=0;
for p=1:npts
    condBit=condBit+sum(data(:,p).^2-2*data(:,p).*expectationsMu.mu-...
        2*(data(:,p)-expectationsMu.mu).*((expectationsB.bChi*...
    expectationsC.c(:,p))))+...
    trace(expectationsB.bbTotalChi*expectationsC.ccT(:,:,p))+...
    trace(expectationsMu.mumuT);
end
priorCBit=trace(expectationsC.ccT(:,:,1));
for i=2:npts
    priorCBit=priorCBit+trace(diag((ones(nTrans,1)-factors.^2).^-1)*...
    (expectationsC.ccT(:,:,i)-...
    2*diag(factors)*expectationsC.cAltc(:,:,i-1)+diag(factors.^2)*expectationsC.ccT(:,:,i-1)));
    
end
priorMuBit=trace(expectationsMu.mumuT);
priorCBit=priorCBit+(npts-1)*sum(log(ones(nTrans,1)-factors.^2));
priorBBit=-nGenes*nTrans*log(alpha)+alpha*trace(expectationsB.bbTotal);
f=nGenes*npts*log(beta)-beta*condBit-priorCBit-priorBBit-...
  priorMuBit+expectationsMu.entropy+expectationsB.entropy+...
  expectationsC.entropy;
