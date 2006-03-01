function f=chipVarLikeGammaGrad(gamma,expectationsC,l)

% CHIPVARLIKEGAMMAGRAD gradient of chipVarLikeGamma

% CHIPVAR
factor=(cos(gamma)+1)/(2+4e-6)+1e-6;
factor=0.99*factor;
npts=size(expectationsC.ccT,3);
preCt=sum(expectationsC.ccT(:,:,2:end),3);
preCtprev=sum(expectationsC.ccT(:,:,1:end-1),3);
preCaltC=sum(expectationsC.cAltc,3);
ct2=preCt(l,l);
ct2prev=preCtprev(l,l);
ctAltct=preCaltC(l,l);
f=(2*factor*((1-factor^2)^-2)*(ct2+factor^2*ct2prev-2*factor*ctAltct)+((1-factor^2)^-1)* ...
   (2*factor*ct2prev-2*ctAltct)-2*(npts-1)*factor/(1-factor^2))*(-0.99*sin(gamma))/(2+4e-6);
