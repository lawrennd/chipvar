% DEMTUVAR demonstrates chipVar on metabolic data

% CHIPVAR
clear all
clear functions
[data,X,annotation,TransNames]=chipVarTuLoadData();
%load dataTu
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
options=chipVarOptions();
%options.stepChecks=1;
model=chipVarInit(X);
randn('seed',198);
pippo=randn(nTrans);
beta=model.beta;
Gamma=model.Gamma;
alpha=model.alpha;
expectationsB.bbTotal=eye(nTrans);
expectationsB.bbTotalChi=eye(nTrans).*(X'*X);
expectationsB.b=randn(nGenes,nTrans);
expectationsB.bChi=expectationsB.b.*X;
expectationsMu.mu=randn(nGenes,1);
expectationsMu.mumuT=eye(size(data,1));
expectationsMu.entropy=0;
expectationsB.entropy=0;%initialises expectations of b
expectationsC=chipVarEstepCMu(data,X,beta,Gamma,expectationsB, expectationsMu);
expectationsB=chipVarEstepBMu(data,X,beta,alpha,expectationsC, expectationsMu);
[model,expectationsB,expectationsC,expectationsMu]= ...
    chipVarEMmu(data,X,model,options,expectationsB,expectationsC, ...
                expectationsMu);

save resultsTu model expectationsB expectationsC expectationsMu