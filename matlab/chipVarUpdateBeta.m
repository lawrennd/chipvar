function beta=chipVarUpdateBeta(data,X,expectationsB,expectationsC,expectationsMu)
% CHIPVARUPDATEBETA updates the precision in E-M algorithm

% CHIPVAR
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
condBit=0;
for p=1:npts
    condBit=condBit+sum(data(:,p).^2-2*data(:,p).*((expectationsB.bChi*...
    expectationsC.c(:,p))))+...
    trace(expectationsB.bbTotalChi*expectationsC.ccT(:,:,p));
end
beta=(condBit/(nGenes*npts))^-1;