CHIPVAR software
Version 0.11		Saturday 06 Jan 2007 at 11:28
Copyright (c) 2007 , 2006, Guido Sanguinetti

Code for inference of transcription factor concentrations and regulatory interactions between genes and transcription factors. 


Version 0.11
------------

Fixed bug where code distribution meant that chipVarEstepCMu.m was an empty file. Thanks to Parantu Shah for bringing this to our attention.

Version 0.1
-----------

This code coincides with the first draft of the paper and recreates the experiments therein.

MATLAB Files
------------

Matlab files associated with the toolbox are:

chipVarEMmu.m: E-M algorithm for dynamical network analysis
chipChipTextRead.m: reads TXT file for the Lee ChIP data files.
chipVarRedLoadData.m: loads reduced Spellman Data with Lee et al ChIP data.
chipVarInit.m: initialises E-M algorithm
chipVarLike.m: likelihood for chipVar to optimise Gamma and m
chipVarUpdateBetaMu.m: updates the precision in E-M algorithm
chipVarEM.m: E-M algorithm for dynamical network analysis
demSpellmanVar.m: demonstrates chipVar on Spellman
chipVarFakeData.m: artificial data for checking chipVar
chipTextRead.m: reads TXT file for the Spellman data files.
chipVarEstepB.m: variational update of the expectations of b
chipVarEstepC.m: variational update of C
chipVarLikeGamma.m: one dimensional gamma likelihood
chipVarEstepBMu.m: variational update of the expectations of b
chipVarLikelihoodBound.m: variational lower bound on marginal likelihood
chipVarEstepCMu.m: variational update of expectations of C.
chipVarLikeGrad.m: gradient of chipVarLike
chipVarLikelihoodCheck.m: Compute the difference in likelhoods.
chipVarTuLoadData.m: loads metabolic Data with combined ChIP data.
chipVarEstep.m: variational E-step for chipVar
chipVarUpdateAlpha.m: updates alpha in E-M
chipTuTextRead.m: assigns common names to probe IDs
chipVarLikelihoodCheckMu.m: Compute the difference in likelhoods.
chipVarUpdateBeta.m: updates the precision in E-M algorithm
chipVarEstepMuMu.m: variational update of mu
chipVarOptions.m: sets default options for chipVar
chipVarLikelihoodBoundMu.m: variational lower bound on marginal likelihood
chipVarEstepMu.m: variational E-step for chipVar
chipVarLikeGammaGrad.m: gradient of chipVarLikeGamma
demTuVar.m: demonstrates chipVar on metabolic data
demFakeVar.m: demonstrates chipVar on artificial data
