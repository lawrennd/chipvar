function f=chipVarLikeGrad(Gamma,data,X,alpha,beta,...
                       expectationsB,expectationsC)
% CHIPVARLIKEGRAD gradient of chipVarLike

% CHIPVAR
nTrans=size(X,2);
npts=size(data,2);
nGenes=size(data,1);
%Gamma=Gammam(1:nTrans)';
%m=Gammam(nTrans+1:end)';
factors=(cos(Gamma)+ones(1,nTrans))/(2+4e-6)+1e-6*ones(1,nTrans);
factors=0.99*factors;
%factors=(1-1e-6)*exp(Gamma)./(ones(1,nTrans)+exp(Gamma));
% Gamma(find(factors==1))=Gamma(find(factors==1))-1e-3;
% factors=cos(Gamma).^2;
preP=zeros(nTrans,npts);
prePQ=zeros(nTrans,npts-1);
for i=1:npts
  preP(:,i)=diag(expectationsC.ccT(:,:,i));
end
for i=1:npts-1
  prePQ(:,i)=diag(expectationsC.cAltc(:,:,i));
end
p=sum(preP(:,2:npts),2)';
q=sum(preP(:,1:npts-1),2)';
pq=sum(prePQ,2)';
gradGamma=+(2*(npts-1)*factors.^3-2*factors.^2.*pq+2*factors.*(p+q- ...
                                                  (npts-1).* ...
                                                  ones(1,nTrans))- ...
           2*pq).*(0.99*sin(Gamma))./((ones(1,nTrans)-factors.^2).^2*(2+4e-6));
%gradOldFactors=(factors-(factors.^2)/(1-1e-6))       
%gradM=(2*(npts-1)*m.*(ones(nTrans,1)-factors).^2-2* ...
%                                    sum(expectationsC.c(:,2:npts), ...
%                                        2).*(ones(nTrans,1)- ...
%                                             factors)+2*factors.* ...
%                                    sum(expectationsC.c(:,1:npts- ...
%                                                  1),2).* ...
%                                    (ones(nTrans,1)-factors))./ ...
%      (ones(nTrans,1)-factors.^2)-2* expectationsC.c(:,1)+2*m;
f=-gradGamma;








