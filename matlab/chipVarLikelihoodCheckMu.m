function [deltaL, L] = chipVarLikelihoodCheckMu(model, expectationsB,...
                                              expectationsC,expectationsMu, data, ...
                                       X, oldL, param);

% CHIPVARLIKELIHOODCHECKMU Compute the difference in likelhoods.

% CHIPVAR
L = chipVarLikelihoodBoundMu(model, data, X,expectationsC,expectationsB,expectationsMu);
deltaL = L-oldL;
oldL = L;
fprintf('Likelihood change of with update of %s: %2.4f\n', param, deltaL)
if deltaL < 0
  warning(['Likelihood drop of ' num2str(deltaL) ' after update of ' param '.']);
end