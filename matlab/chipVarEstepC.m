function expectationsC=chipVarEstepC(data,X,beta,Gamma, expectationsB)
% CHIPVARESTEPC variational update of C

% CHIPVAR
nTrans=size(X,2);
nGenes=size(data,1);
npts=size(data,2);
factors=(cos(Gamma)+ones(nTrans,1))/(2+4e-6)+1e-6*ones(nTrans,1);
factors=0.99*factors;
%factors=(1-1e-6)*exp(Gamma)./(ones(nTrans,1)+exp(Gamma));
% Gamma(find(factors==1))=Gamma(find(factors==1))-1e-3;
% factors=cos(Gamma).^2;
expectationsB.bChi=expectationsB.b.*X;
block1=diag((ones(nTrans,1)-factors.^2).^-1);
block2=diag(ones(nTrans,1)+factors.^2);
K=[block1+beta*expectationsB.bbTotalChi,-diag(factors)*block1,...
    zeros(nTrans,(npts-2)*nTrans)];
for i=2:npts-1
    K=[K;zeros(nTrans,(i-2)*nTrans),-diag(factors)*block1,block2*block1+...
        beta*expectationsB.bbTotalChi,-diag(factors)*block1,...
    zeros(nTrans,(npts-1-i)*nTrans)];
end
K=[K;zeros(nTrans,(npts-2)*nTrans),-diag(factors)*block1,...
    block1+beta*expectationsB.bbTotalChi];
expectationsC.entropy=logdet(K);
invDiagCholBlock=zeros(nTrans,nTrans,npts);
superDiagCholBlock=zeros(nTrans,nTrans,npts-1);
invDiagCholBlock(:,:,1)=inv(chol(block1+beta*expectationsB.bbTotalChi));
for i=2:npts-1
  superDiagCholBlock(:,:,i-1)=invDiagCholBlock(:,:,i-1)'*(- ...
                                                    diag(factors)*block1);
  invDiagCholBlock(:,:,i)=inv(chol(block2*block1+...
        beta*expectationsB.bbTotalChi-superDiagCholBlock(:,:,i-1)'* ...
                                     superDiagCholBlock(:,:,i-1)));
end
superDiagCholBlock(:,:,npts-1)=invDiagCholBlock(:,:,npts-1)'*(- ...
                                                    diag(factors)*block1);
invDiagCholBlock(:,:,npts)=inv(chol(block1+...
        beta*expectationsB.bbTotalChi-superDiagCholBlock(:,:,npts-1)'* ...
                                     superDiagCholBlock(:,:,npts- ...
                                                  1)));
postcov.diag=zeros(nTrans,nTrans,npts);
postcov.superDiag=zeros(nTrans,nTrans,npts-1);
postcov.diag(:,:,npts)=invDiagCholBlock(:,:,end)*invDiagCholBlock(:,:,end)';
for i=1:npts-1
  postcov.superDiag(:,:,npts-i)=-invDiagCholBlock(:,:,npts-i)* ...
      superDiagCholBlock(:,:,npts-i)*postcov.diag(:,:,npts-i+1);
  postcov.diag(:,:,npts-i)=invDiagCholBlock(:,:,npts-i)* ...
      invDiagCholBlock(:,:,npts-i)'-postcov.superDiag(:,:,npts-i)* ...
      superDiagCholBlock(:,:,npts-i)'*invDiagCholBlock(:,:,npts-i)';
end
postCov=zeros(nTrans*npts,nTrans*npts);
for i=1:npts
  postCov((i-1)*nTrans+1:i*nTrans,(i-1)*nTrans+1:i*nTrans)= ...
      postcov.diag(:,:,i)/2;
end
for i=1:npts-1
  postCov((i-1)*nTrans+1:i*nTrans,i*nTrans+1:(i+1)*nTrans)= ...
      postcov.superDiag(:,:,i);
end
for i=1:npts-2
  for j=npts-i+1:npts
    postCov((npts-i-2)*nTrans+1:(npts-i-1)*nTrans,(j-1)*nTrans+1:j* ...
            nTrans)=postCov((npts-i-2)*nTrans+1:(npts-i-1)*nTrans, ...
                            (npts-i-1)*nTrans+1:(npts-i)*nTrans)* ...
        pdinv(postCov((npts-i-1)*nTrans+1:(npts-i)*nTrans,(npts- ...
                                                      i-1)*nTrans+ ...
                      1:(npts-i)*nTrans))*postCov((npts-i-1)* ...
                                                    nTrans+1:(npts- ...
                                                      i)*nTrans, ...
                                                    (j-1)*nTrans+1:j*nTrans)/2;
  end
end
postCov=(postCov+postCov');
bBit=beta*((data)'*expectationsB.bChi)';

bBit=bBit(:);
%preExpectations=triblocksolve(K,bBit,npts);
preExpectations=postCov*bBit;
expectationsC.c=reshape(preExpectations,nTrans,npts);

expectationsC.ccT=zeros(nTrans,nTrans,npts);
for t=1:npts
    expectationsC.ccT(:,:,t)=postcov.diag(:,:,t)+expectationsC.c(:,t)*expectationsC.c(:,t)';
end
expectationsC.cAltc=zeros(nTrans,nTrans,npts-1);
for t=1:npts-1
    expectationsC.cAltc(:,:,t)=postcov.superDiag(:,:,t)+ ...
        expectationsC.c(:,t)*expectationsC.c(:,t+1)';
    
end
