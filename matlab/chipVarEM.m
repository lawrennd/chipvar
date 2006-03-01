function [model,expectationsB,expectationsC]=...
    chipVarEM(data,X,model,options,expectationsB,expectationsC)
% CHIPVAREM E-M algorithm for dynamical network analysis

% CHIPVAR

%randn('seed',49);
nGenes=size(X,1);
nTrans=size(X,2);
npts=size(data,2);
params = {'alpha', 'beta', 'Gamma', 'estep'};
oldL=chipVarLikelihoodBound(model,data,X,...
                      expectationsC, expectationsB);
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
      model.beta = chipVarUpdateBeta(data, X, expectationsB, expectationsC);
     
     case 'Gamma'
      options.optimizer(14)=100;
      %options.optimizer(9)=1;
      %options.optimizer(1)=1;
      Gamma=model.Gamma';
      Gamma=scg('chipVarLike',Gamma,options.optimizer,'chipVarLikeGrad', ...
                 data,X,model.alpha,model.beta,expectationsB,expectationsC);
      
      model.Gamma=Gamma';
%     case 'm'
%      model.m=chipVarUpdateM(model,expectationsC);
     case 'estep'
      [expectationsB,expectationsC] = chipVarEstep(data,X,model,...
                                                   expectationsB,expectationsC,2,options);  
    end
    if options.stepChecks
      [deltaL, oldL] = chipVarLikelihoodCheck(model, expectationsB,...
                                              expectationsC, data,X, oldL,...
                                              updateParam);
      maxDeltaL = max([maxDeltaL deltaL]);
    end
  end
  if ~options.stepChecks
    [maxDeltaL, oldL] = chipVarLikelihoodCheck(model, expectationsB,...
                                               expectationsC, data,X, oldL,...
                                               'All params');
  end
  lastUpdate = order(end);
  order = randperm(length(params));
  while(lastUpdate == order(1))
    order = randperm(length(params));
  end
   if rem(counter,10)==0
    save partialResultsSpellmanFast
  end
  fprintf('Iteration number: %d\n', counter);
end
