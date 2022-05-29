function mutatedOffspring = mutation(toBeMutated)

%% Mutation function
% This function performs the inversion of a random allele in the offspring
% genome with a default probability of mutProb.

global mutProb

genomeLength = length(toBeMutated); % Number of alleles in the genome
mutatedAllele = randperm(genomeLength,1); % Random selection of the allele to be inverted
    
% Mutation
    if rand <= mutProb                   
        mutatedOffspring = toBeMutated;        
        if toBeMutated(mutatedAllele) == 1
            mutatedOffspring(mutatedAllele) = 0;            
        else
            mutatedOffspring(mutatedAllele) = 1;
        end          
        
    else                                        
        mutatedOffspring = toBeMutated;    
    end
end