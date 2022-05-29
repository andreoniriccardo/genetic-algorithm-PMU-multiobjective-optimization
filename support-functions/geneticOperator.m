function childOffspring  = geneticOperator(parentsSelected)
global G
%% Genetic operator
% The crossover process is performed
% The input is the matrix of the selected parents
% The output is the matrix of the offsprings
%% Crossover

[parentSelSize] = size(parentsSelected,1);
rc = randperm(parentSelSize); % Random permutation
childOffspring = zeros(parentSelSize,G);

for i = 1:(parentSelSize/2)
    % Parent couples are randomly selected
    parent1 = parentsSelected((rc(2*i-1)),:);
    parent2 = parentsSelected((rc(2*i)),:);
    % If parents are equal, so are offsprings
    if (isequal(parent1,parent2)) == 1 && rand(1) > 0.5

        child1 = parent1;
        child2 = parent2;
    % Otherwise binary TPC is performed
    else

        temp = 1:1:G+1;
        randIdcs = sort(randperm(length(temp),4));
        child1 = [parent1(1:(randIdcs(1)-1)) parent2(randIdcs(1):(randIdcs(2)-1)) parent1(randIdcs(2):(randIdcs(3)-1)) parent2(randIdcs(3):(randIdcs(4)-1)) parent1(randIdcs(4):G)];
        child2 = [parent2(1:(randIdcs(1)-1)) parent1(randIdcs(1):(randIdcs(2)-1)) parent2(randIdcs(2):(randIdcs(3)-1)) parent1(randIdcs(3):(randIdcs(4)-1)) parent2(randIdcs(4):G)];

    end

    % Offspring matrix
    childOffspring((rc(2*i-1)),:) = mutation(child1);
    childOffspring((rc(2*i)),:) = mutation(child2);
end
