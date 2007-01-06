function [geneName, annotation, data] = chipChipTextRead(file1, file2)

% CHIPCHIPTEXTREAD reads TXT file for the Lee ChIP data files.

% CHIPVAR
[geneName,annotation]=...
    textread(file1,'%q %q %*[^\n]',...
    'headerlines',2,'whitespace','','delimiter','\t');

data = load(file2);