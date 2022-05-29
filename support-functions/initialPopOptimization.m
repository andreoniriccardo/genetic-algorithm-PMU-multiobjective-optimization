function [initialPop] = initialPopOptimization()

global N G Incidence

% Trovo la soluzione con il minimo numero di PMU
solMinPMU = intlinprog(ones(1,G),G,-Incidence,-ones(G,1),[],[],zeros(G,1),ones(G,1));
% Trovo il minimo numero di PMU
minPMU = nnz(solMinPMU);

totPart = round(N/(G-minPMU+1)-1);

intcon = G;
b = ones(G,1);
f = zeros(1,G);
Aeq = ones(1,G);
lb = zeros(G,1);
ub = ones(G,1);

count = 1; % Tiene il conto del totale delle soluzioni generate + 1
for i=minPMU:G-1 % Per qualunque possibile valore di n PMU ammissibile
    if i == minPMU
        pool(count,:) = solMinPMU';
        pool = round(pool);
        disp({'Trovata soluzione con _ PMU',i});
        disp({'Totale soluzioni trovate:',count});
        temp = pool(count,:);
        count = count + 1;
    else
        % Determino una possibile soluzione
        pool(count,:) = intlinprog(f,intcon,-Incidence,-b,Aeq,i,lb,ub);
        pool = round(pool);
        disp({'Trovata soluzione con _ PMU',i});
        disp({'Totale soluzioni trovate:',count});
        temp = pool(count,:);
        count = count + 1;
    end
    for j = 1:(totPart-1) % Determino totPart-1 (una � gi� stata appena creata) possibili soluzioni ammissibili
        % Seleziono casualmente un cromosoma 1 dall'ultima soluzione generata
        clear oneIndices
        if nnz(temp) < G
            oneIndices = find(~temp); % Trovo gli indici degli alleli pari a 0
            randomIndex = randsample(oneIndices,1); % Seleziono casualmente uno di questi indici
        else
            randomIndex = randsample(1:G,1);
        end

        % Vario i vincoli di uguaglianza per imporre che tale indice sia
        % pari a 0
        clear newConstr
        newConstr = zeros(1,G);
        newConstr(randomIndex) = 1;
        Aeq_new = [ Aeq;
                    newConstr];
        beq_new = [i;
                   1];
        % Calcolo la nuova soluzione con i PMU
        soluzione = intlinprog(f,intcon,-Incidence,-b,Aeq_new,beq_new,lb,ub);
        soluzione = soluzione';
        solSize = size(soluzione);
        if solSize(2) == G
            pool(count,:) = soluzione;
            pool = round(pool);
            disp({'Trovata soluzione con _ PMU (+)',i});
            disp({'Totale soluzioni trovate:',count});
            temp = pool(count,:);
            count = count + 1;
        else
            pool(count,:) = ones(1,G);
            pool = round(pool);
            disp({'Trovata soluzione con _ PMU (-)',i});
            disp({'Totale soluzioni trovate:',count});
            temp = pool(count,:);
            count = count + 1;
        end
    end

end

a = size(pool);
pmu30 = ones(N-a(1),G);
pool=[  pool;
        pmu30];

initialPop = round(pool);


end
