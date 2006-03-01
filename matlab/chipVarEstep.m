function [expectationsB,expectationsC]=chipVarEstep(data,X,model, expectationsB,expectationsC,maxCounter,options)
% CHIPVARESTEP variational E-step for chipVar

% CHIPVAR

nGenes=size(data,1);
npts=size(data,2);
beta=model.beta;
alpha=model.alpha;
Gamma=model.Gamma;
change=1;
counter=0;
initialLike=chipVarLikelihoodBound(model,data,X,...
                      expectationsC, expectationsB);
while counter < maxCounter
    counter=counter+1;
    expectationsC=chipVarEstepC(data,X,beta,Gamma,expectationsB);
    if options.stepChecks
       likeBound=chipVarLikelihoodBound(model,data,X,...
                      expectationsC, expectationsB);
      
    if counter>1
          change1=likeBound-newLikeBound;
           fprintf('Likelihood change of %d with update of c\n',change1)
          %change=max(change1,change2);
          if change1<0
              warning('Likelihood drop of %d after update of C\n',change1)
              counter
          end
    else 
      change1=likeBound-initialLike;
       fprintf('Likelihood change of %d with update of c\n',change1)
      if change1<0
              warning('Likelihood drop of %d after update of c\n',change1)
              counter
          end
    end
    end
    expectationsB=chipVarEstepB(data,X,beta,alpha, expectationsC);
    if options.stepChecks
    newLikeBound=chipVarLikelihoodBound(model,data,X,...
                      expectationsC, expectationsB); 
    change2=newLikeBound-likeBound;
     fprintf('Likelihood change of %d with update of b\n',change2)
    if change2<-1e-3
          warning('Likelihood drop of %d after update of b\n',change2)
          counter
    end
    end
end
%finalLike=max([newLikeBound,likeBound]);
