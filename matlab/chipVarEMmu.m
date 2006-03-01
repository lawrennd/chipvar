function [model,expectationsB,expectationsC,expectationsMu]=...
    chipVarEMmu(data,X,model,options,expectationsB,expectationsC,expectationsMu)
% CHIPVAREMMU E-M algorithm for dynamical network analysis

% CHIPVAR
rand('seed',90);
nGenes=size(X,1);
nTrans=size(X,2);
npts=size(data,2);
beta=model.beta;
Gamma=model.Gamma;
alpha=model.alpha;
params = {'alpha', 'beta', 'Gamma', 'estep'};
oldL=chipVarLikelihoodBoundMu(model,data,X,...
                      expectationsC, expectationsB,expectationsMu);
maxDeltaL=1;
counter=0;
while (maxDeltaL > options.tol & counter < options.maxIters)
  maxDeltaL = 0;
  counter=counter+1;
  if counter == 1
    order = randperm(length(params));
  end
  for i = order
    updateParam = params{i};
    switch updateParam
     case 'alpha'
%      model.mu = chipVarUpdateMu(expectationsB);
      model.alpha = chipVarUpdateAlpha(expectationsB);
      
     case 'beta'
      model.beta = chipVarUpdateBetaMu(data, X, expectationsB, expectationsC,expectationsMu);
     
     case 'Gamma'
      options.optimizer(14)=10;
      %options.optimizer(9)=1;
      %options.optimizer(1)=1;
      Gamma=model.Gamma';
      for l=1:nTrans 
        Gamma(l)=scg('chipVarLikeGamma',Gamma(l), ...
                     options.optimizer,'chipVarLikeGammaGrad', expectationsC,l);
        
      end
      model.Gamma=Gamma';
%     case 'm'
%      model.m=chipVarUpdateM(model,expectationsC);
     case 'estep'
      [expectationsB,expectationsC,expectationsMu] = ...
          chipVarEstepMu(data,X,model,expectationsB,expectationsC,expectationsMu,2,options);
      
    end
    if options.stepChecks
      [deltaL, oldL] = chipVarLikelihoodCheckMu(model, expectationsB,...
                                              expectationsC,expectationsMu, data,X, oldL,...
                                              updateParam);
      maxDeltaL = max([maxDeltaL deltaL]);
    end
  end
  if ~options.stepChecks
    [maxDeltaL, oldL] = chipVarLikelihoodCheckMu(model, expectationsB,...
                                               expectationsC, expectationsMu,data,X, oldL,...
                                               'All params');
  end
  lastUpdate = order(end);
  order = randperm(length(params));
  while(lastUpdate == order(1))
    order = randperm(length(params));
  end
  fprintf('Iteration number: %d\n', counter);
  if rem(counter,10)==0
    save partialResultsTu
  end
end
