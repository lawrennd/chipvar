CHIPVAR software
Version 0.1		Monday 22 May 2006 at 21:51
Copyright (c) 2006 Neil D. Lawrence

Code for inference of transcription factor concentrations and regulatory interactions between genes and transcription factors. 

Version 0.1
-----------

This code coincides with the first draft of the paper and recreates the experiments therein.

MATLAB Files
------------

Matlab files associated with the toolbox are:

chipChipTextRead.m: reads TXT file for the Lee ChIP data files.
chipTextRead.m: reads TXT file for the Spellman data files.
chipTuTextRead.m: assigns common names to probe IDs
chipVarEM.m: E-M algorithm for dynamical network analysis
chipVarEMmu.m: E-M algorithm for dynamical network analysis
chipVarEstepB.m: variational update of the expectations of b
chipVarEstepBMu.m: variational update of the expectations of b
chipVarEstepC.m: variational update of C
chipVarEstep.m: variational E-step for chipVar
chipVarEstepMu.m: variational E-step for chipVar
chipVarEstepMuMu.m: variational update of mu
chipVarFakeData.m: artificial data for checking chipVar
chipVarInit.m: initialises E-M algorithm
demTuVar.m: demonstrates chipVar on metabolic data
chipVarLikeGammaGrad.m: gradient of chipVarLikeGamma
chipVarLikeGamma.m: one dimensional gamma likelihood
chipVarLikeGrad.m: gradient of chipVarLike
chipVarLikelihoodBound.m: variational lower bound on marginal likelihood
chipVarLikelihoodBoundMu.m: variational lower bound on marginal likelihood
chipVarLikelihoodCheck.m: Compute the difference in likelhoods.
chipVarLikelihoodCheckMu.m: Compute the difference in likelhoods.
chipVarLike.m: likelihood for chipVar to optimise Gamma and m
chipVarOptions.m: sets default options for chipVar
chipVarRedLoadData.m: loads reduced Spellman Data with Lee et al ChIP data.
chipVarTuLoadData.m: loads metabolic Data with combined ChIP data.
chipVarUpdateAlpha.m: updates alpha in E-M
chipVarUpdateBeta.m: updates the precision in E-M algorithm
chipVarUpdateBetaMu.m: updates the precision in E-M algorithm
demFakeVar.m: demonstrates chipVar on artificial data
demSpellmanVar.m: demonstrates chipVar on Spellman
