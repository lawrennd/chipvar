function [expectationsB,expectationsC,expectationsMu]= ...
    chipVarEstepMu(data,X,model, expectationsB,expectationsC,expectationsMu,maxCounter,options)

% CHIPVARESTEPMU variational E-step for chipVar

% CHIPVAR
nGenes=size(data,1);
npts=size(data,2);
beta=model.beta;
alpha=model.alpha;
Gamma=model.Gamma;
change=1;
counter=0;
while counter < maxCounter
    counter=counter+1;
  %  warning('pippo')
    expectationsMu=chipVarEstepMuMu(data,beta,expectationsB, ...
                                  expectationsC);
    if options.stepChecks
      likeBound1=chipVarLikelihoodBoundMu(model,data,X,...
                      expectationsC, expectationsB,expectationsMu);
      if counter >1 
      change1=likeBound1-likeBound3;
          %change=max(change1,change2);
          fprintf('Likelihood change of %d with update of Mu\n',change1)
          if change1<0
              warning('Likelihood drop of %d after update of mu\n',change1)
              counter
          end
      end
      
    end
    expectationsB=chipVarEstepBMu(data,X,beta,alpha, expectationsC, ...
                                  expectationsMu);
    if options.stepChecks
    likeBound2=chipVarLikelihoodBoundMu(model,data,X,...
                      expectationsC, expectationsB,expectationsMu); 
    change2=likeBound2-likeBound1;
    fprintf('Likelihood change of %d with update of b\n',change2)
    if change2<0
          warning('Likelihood drop of %d after update of b\n',change2)
          counter
    end
    end
   
   expectationsC=chipVarEstepCMu(data,X,beta,Gamma,expectationsB, ...
                                 expectationsMu);
   if options.stepChecks
   likeBound3=chipVarLikelihoodBoundMu(model,data,X,...
                     expectationsC, expectationsB,expectationsMu);
         change3=likeBound3-likeBound2;
         fprintf('Likelihood change of %d with update of c\n',change3)
         %change=max(change1,change2);
         if change3<0
             warning('Likelihood drop of %d after update of C\n',change3)
             counter
         end
   end
    
end

