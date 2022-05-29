function [parentSelected] = tournamentSelection(population)
%% Tournament Selection

[popSize, idxCD] = size(population);
idxFront = idxCD - 1;          

% Random permutation of parent population
tournament = [randperm(popSize); randperm(popSize)]';

% Initialization of output matrix
parentSelected = zeros(popSize, idxCD);

% Selecting parents according to the non-dominated-CD criteria
for i = 1: popSize
    couple = tournament(i,:);
 if population(couple(1),idxFront) ~= population(couple(2),idxFront)  
    if population(couple(1),idxFront) < population(couple(2),idxFront)                                                                                                  
        candidateMin = population(couple(1),:);                 
    elseif population(couple(1),idxFront) > population(couple(2),idxFront) 
        candidateMin = population(couple(2),:);     
    end
 parentSelected(i,:) = candidateMin;  
 else           
    if population(couple(1),idxCD ) > population(couple(2),idxCD)
        candidateMax = population(couple(1),:);   
    elseif population(couple(1),idxCD) < population(couple(2),idxCD)
        candidateMax = population(couple(2),:);   
    else
        temp = randperm(2);
        candidateMax = population(couple(temp(1)),:);
    end 
 parentSelected(i,:) = candidateMax; 
 end
end
