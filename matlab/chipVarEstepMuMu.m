function expectationsMu=chipVarEstepMuMu(data,beta,expectationsB, ...
                                         expectationsC)
% CHIPVARESTEPMUMU variational update of mu

% CHIPVAR
npts=size(data,2);
nGenes=size(data,1);
%expectationsB.bChi=expectationsB.b.*X;
postCov=(1+npts*beta)^-1*eye(nGenes);
expectationsMu.mu=beta*((1+npts*beta)^-1)*(sum(data,2)- ...
                                         sum(expectationsB.bChi* ...
                                             expectationsC.c,2));
expectationsMu.mumuT=postCov+expectationsMu.mu*expectationsMu.mu';
expectationsMu.entropy=-nGenes*log(1+npts*beta);
