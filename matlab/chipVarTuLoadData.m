function [data,precs,X,annotation,TransNames]=chipVarTuLoadData();
% CHIPVARTULOADDATA loads metabolic Data with combined ChIP data.

% CHIPVAR
[probeName, data, vars] = chipTuTextRead;
data=data(find(sum(vars,2)),:);
probeName=probeName(find(sum(vars,2)));
vars=vars(find(sum(vars,2)),:);
[probeName1, annota1, dataChip1] = chipChipTextRead(['../data/' ...
                    'Yeast_Connectivity.txt'], '../data/Connectivity_Matrix.txt');
preDataChip=load('../data/Connectivity2.txt');
[preProbeName2, preAnnota] = textread('../data/annotations2.txt','%q%q');
TransNames=textread('../data/Trans_Names2.txt','%q');
TransNames1=textread('../data/Trans_Names.txt','%q', ...
                    'headerlines',1,'whitespace','','delimiter','\t');
TransNames1=TransNames1(4:end-2);
%TransNames=TransNames(2:end);
counter=0;
dataChip=[];
for i=1:size(dataChip1,1)
       gigio=strcmp(probeName1(i),preProbeName2);
       if sum(gigio)==1
         dataChip=[dataChip;preDataChip(find(gigio),:)];
       else
         counter=counter+1;
         dataChip=[dataChip;ones(1,size(TransNames,1))];
       end
end
dataChip2=[];
for i=1:size(TransNames,1)
  pippo=strcmp(TransNames(i),TransNames1);
  if sum(pippo)==1
    dataChip2=[dataChip2,dataChip1(:,find(pippo))];
  else
    dataChip2=[dataChip2,ones(size(dataChip1,1),1)];
  end
end



index=zeros(size(dataChip,1),1);
for i=1:size(dataChip,1)
    index(i)=sum(strcmp(probeName1(i),probeName));
   
end
dataChip=dataChip(find(index),:);
dataChip2=dataChip2(find(index),:);
annota1=annota1(find(index));
probeName1=probeName1(find(index));
index=zeros(size(data,1),1);
preX1=[];
preX2=[];
annotation=[];
for i=1:size(data,1)
    index(i)=sum(strcmp(probeName(i),probeName1));
    if index(i)
      preX1=[preX1;dataChip(find(strcmp(probeName(i), probeName1)), ...
                            :)];
      preX2=[preX2;dataChip2(find(strcmp(probeName(i), probeName1)),:)];
      annotation=[annotation; annota1(find(strcmp(probeName(i),probeName1)))];
    end
end
data=data(find(index),:);
vars=vars(find(index),:);
probeName=probeName(find(index));
X1=zeros(size(preX1,1),size(preX1,2));
X2=zeros(size(preX1,1),size(preX1,2));
I1=find(preX1<1e-3);
X1(I1)=1;
I2=find(preX2<1e-3);
X2(I2)=1;
X=zeros(size(X1));
X(find(X1+X2))=1;
%X=X(:,1:20);

fakeX=sum(X,2);
X=X(find(fakeX),:);
annotation=annotation(find(fakeX));
effectX=sum(X,1);
TransNames=TransNames(find(effectX));
X=X(:,find(effectX));
data=data(find(fakeX),:);
%vars=vars(find(fakeX),:);
%precs=vars.^-2;