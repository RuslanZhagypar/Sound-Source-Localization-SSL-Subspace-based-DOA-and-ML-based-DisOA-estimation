clear all
% outputs = [[2,1],[2,1],[3,1]];
% targets = [[4,0],[5,0],[4,0]];
outputs = [0 1 1];
targets = [1 0 1];
[c,cm] = confusion(outputs,targets);
fprintf('Percentage Correct Classification   : %f%%\n', 100*(1-c));