function [deltaL, L] = chipVarLikelihoodCheck(model, expectationsB,...
                                              expectationsC, data, ...
                                       X, oldL, param);

% CHIPVARLIKELIHOODCHECK Compute the difference in likelhoods.

% CHIPVAR
L = chipVarLikelihoodBound(model, data, X,expectationsC,expectationsB);
deltaL = L-oldL;
oldL = L;
fprintf('Likelihood change of with update of %s: %2.4f\n', param, deltaL)
if deltaL < 0
  warning(['Likelihood drop of ' num2str(deltaL) ' after update of ' param '.']);
end