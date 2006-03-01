function expectationsB=chipVarEstepB(data,X,beta,alpha,expectationsC)
% CHIPVARESTEPB variational update of the expectations of b

% CHIPVAR
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
fixedBit=beta*sum(expectationsC.ccT,3)/alpha;
expectationsB.b=[];
postCov=zeros(nTrans,nTrans);
postCovChi=zeros(nTrans,nTrans);
expectationsB.entropy=0;
for i=1:nGenes
    if alpha>0.1
      numX=sum(X(i,:));
      YMat=zeros(nTrans,numX);
      [crap, index]=find(X(i,:));
      for j=1:numX
         YMat(index(j),j)=1;
      end


    
      auxMat=fixedBit-fixedBit*YMat*pdinv(eye(size(YMat,2))+YMat'* ...
                                              fixedBit*YMat)* ...
       YMat'*fixedBit;
      newBit=alpha^-1*(eye(nTrans)-auxMat.*(X(i,:)'*X(i,:)));
    else
      newBit=pdinv(alpha*eye(nTrans)+fixedBit.*(X(i,:)'*X(i,:)));
    end
    
    postCov=postCov+newBit;
    expectationsB.entropy=expectationsB.entropy-logdet(newBit);
    newBitChi=newBit.*(X(i,:)'*X(i,:));
    postCovChi=postCovChi+newBitChi;
    expectationsB.b=[expectationsB.b;beta*((newBit* ...
                                              (expectationsC.c* ...
                                               data(i,:)')).*X(i,:)')'];
end
expectationsB.bbTotal=postCov+expectationsB.b'*expectationsB.b;
expectationsB.bbTotalChi=postCovChi+(expectationsB.b.*X)'*...
    (expectationsB.b.*X);
expectationsB.bChi=expectationsB.b.*X;    