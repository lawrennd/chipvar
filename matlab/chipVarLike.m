function f=chipVarLike(Gamma,data,X,alpha,beta,...
                       expectationsB,expectationsC)
% CHIPVARLIKE likelihood for chipVar to optimise Gamma and m

% CHIPVAR
nTrans=size(X,2);
model.Gamma=Gamma';
model.alpha=alpha;
model.beta=beta;
f=-chipVarLikelihoodBound(model,data,X,expectationsC, expectationsB);