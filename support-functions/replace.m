function newPop = replace(interPopulationSorted)
global N G F
%% Replace function
% Select N individuals among the intermediate population (of dimension 2*N)
% Select individuals belonging to less dominated fronts
% if the selection of an entire front exceeds the total of N individuals,
% select the individuals in that front that have higher values of CD

temp = sortrows(interPopulationSorted, [G+F+2, G+F+3],{'ascend' 'descend'});
newPop = temp(1:N,:);

end