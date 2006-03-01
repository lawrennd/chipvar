
% DEMFAKEVAR demonstrates chipVar on artificial data

% CHIPVAR

clear all
clear functions
[data,X,annotation,TransNames]=chipVarRedLoadData();
[data,originalC, originalB, originalGamma, originalAlpha,originalBeta]=chipVarFakeData(X,8);
nGenes=size(data,1);
npts=size(data,2);
nTrans=size(X,2);
options=chipVarOptions();
%options.stepChecks=1;
model=chipVarInit(X);
beta=model.beta;
Gamma=model.Gamma;
alpha=model.alpha;
%pippo=randn(nTrans);
expectationsB.bbTotal=eye(nTrans);
expectationsB.bbTotalChi=eye(nTrans).*(X'*X);
expectationsB.b=randn(nGenes,nTrans);
expectationsB.bChi=expectationsB.b.*X;
expectationsB.entropy=0;%initialises expectations of b
expectationsC=chipVarEstepC(data,X,beta,Gamma,expectationsB);
expectationsB=chipVarEstepB(data,X,beta,alpha,expectationsC);
[model,expectationsB,expectationsC]=chipVarEM(data,X,model,options, ...
                                              expectationsB,expectationsC);
save resultsFake model expectationsB expectationsC
