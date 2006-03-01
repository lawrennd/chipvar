function [ORF, data, vars]=chipTuTextRead()
% CHIPTUTEXTREAD assigns common names to probe IDs

% CHIPVAR
[IDs,crap,ORF,geneCN]=textread('../data/MetabolData/dictionary.txt','%q%q%q%q');
IDSnew=textread('../data/MetabolData/probeIDTu.txt','%q');
data=load('../data/MetabolData/YeastMetabolism_exprs.txt');
vars=load('../data/MetabolData/YeastMetabolism_se.txt');
