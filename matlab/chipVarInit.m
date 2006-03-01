function model=chipVarInit(X)
% CHIPVARINIT initialises E-M algorithm

% CHIPVAR
randn('seed',38);
nGenes=size(X,1);
nTrans=size(X,2);
model.beta=30;
model.alpha=3;
model.Gamma=0.1*randn(nTrans,1);
