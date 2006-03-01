function options=chipVarOptions()
% CHIPVAROPTIONS sets default options for chipVar

% CHIPVAR
options.tol=1e-4;
options.maxIters=1500;
options.tolEstep=1e-4;
options.stepChecks=0;
options.optimizer=foptions;