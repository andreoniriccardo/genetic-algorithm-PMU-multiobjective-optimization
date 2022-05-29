function [initialPop] = initialPopOpt()

global N G Incidence

minPMU = 10;
totPart = round(N/(G-minPMU+1)-1);

intcon = G;
b = ones(G,1);
f = zeros(1,G);
Aeq = ones(1,G);
lb = zeros(G,1);
ub = ones(G,1);

count = 1; % Tiene il conto del totale delle soluzioni generate + 1
for i=minPMU:G-1 % Per qualunque possibile valore di n PMU ammissibile
    
    % Determino una possibile soluzione
    pool(count,:) = intlinprog(f,intcon,-Incidence,-b,Aeq,i,lb,ub);
    pool = round(pool);
    temp = pool(count,:);
    count = count + 1;
    for j = 1:(totPart-1) % Determino totPart-1 (una è già stata appena creata) possibili soluzioni ammissibili
        % Seleziono casualmente un cromosoma 1 dall'ultima soluzione generata
        clear oneIndices
        oneIndices = find(temp); % Trovo gli indici degli alleli pari a 1
        randomIndex = randsample(oneIndices,1); % Seleziono casualmente uno di questi indici
        
        % Vario i vincoli di uguaglianza per imporre che tale indice sia
        % pari a 0
        clear newConstr
        newConstr = zeros(1,G);
        newConstr(randomIndex) = 1;
        Aeq_new = [ Aeq;
                    newConstr];
        beq_new = [i;
                   0];
        % Calcolo la nuova soluzione con i PMU
        pool(count,:) = intlinprog(f,intcon,-Incidence,-b,Aeq_new,beq_new,lb,ub);
        pool = round(pool);
        temp = pool(count,:);
        count = count + 1;
    end
    
end

a = size(pool);
pmu30 = ones(N-a(1),G);
pool=[  pool;
        pmu30];

initialPop = round(pool);

%{
% Per verifica
aa=size(pool);
for iii = 1:aa(1)
    pool(iii,G+1) = nnz(pool(iii,:));
end
%}

end