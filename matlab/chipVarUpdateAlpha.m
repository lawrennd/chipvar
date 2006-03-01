function alpha=chipVarUpdateAlpha(expectationsB)
% CHIPVARUPDATEALPHA updates alpha in E-M

% CHIPVAR
nGenes=size(expectationsB.b,1);
nTrans=size(expectationsB.b,2);
alpha=nGenes*nTrans/trace(expectationsB.bbTotal);