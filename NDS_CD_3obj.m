function [populationSorted] = NDS_CD_3obj(population)
%% Non-dominated sorting and crowding distance computation of a 3 objective problem
    global G F problemType

    %% Iniziallizzazione delle variabili
    
    dominationSet.sp = [];
    frontCount = 1;
    populationSorted1 = [];
    infeasiblePop = [];
    front.fr = [];
    
    %% Feasible solutions identification
    % 3 possible cases
    if all(population(:,G+F+1)==0)        % All solution feasible 
        problemType = 'f';                   
        feasSolutions = population(:,1:G+F);       
        feasSize = size(feasSolutions,1);      
    elseif all(population(:,G+F+1)~=0)    % All solution unfeasible   
        problemType = 'u';                  
        feasSize = 0;                     
        infSolutions = population;    % Mixed  
    else
        problemType = 'm';
        feasIndex = find(population(:,G+F+1)==0);     
        feasSolutions = population(feasIndex,1:G+F);    
        feasSize = size(feasSolutions,1);             
        infeasIndex = find(population(:,G+F+1)~=0);   
        infSolutions = population(infeasIndex,1:G+F+1); 
    end
    
    
    %%
    if problemType=='f' || problemType=='m'   
        
        ObjF1 = feasSolutions(:,G+1);        
        ObjF2 = feasSolutions(:,G+2);       
        ObjF3 = feasSolutions(:,G+3);        

        %% Non-Dominating Sorting 
        for p=1:feasSize
            % A dominates B if: 
            % f1(A)<f1(B) AND f2(A)<f2(B) AND f3(A)<f3(B) |OR| 
            % f1(A)=f1(B) AND f2(A)<f2(B) AND f3(A)<f3(B) |OR|
            % f1(A)=f1(B) AND f2(A)=f2(B) AND f3(A)<f3(B) |OR|
            % f1(A)<f1(B) AND f2(A)=f2(B) AND f3(A)<f3(B) |OR|
            % f1(A)<f1(B) AND f2(A)=f2(B) AND f3(A)=f3(B) |OR|
            % f1(A)<f1(B) AND f2(A)<f2(B) AND f3(A)=f3(B) |OR|
            % f1(A)=f1(B) AND f2(A)<f2(B) AND f3(A)=f3(B)
            
            dominationSet(p).sp = find(((ObjF1(p)-ObjF1)<0  & (ObjF2(p)-ObjF2)<0  & (ObjF3(p)-ObjF3)<0 )  | ...
                                ((ObjF1(p)-ObjF1)==0 & (ObjF2(p)-ObjF2)<0  & (ObjF3(p)-ObjF3)<0 )  | ...
                                ((ObjF1(p)-ObjF1)==0 & (ObjF2(p)-ObjF2)==0 & (ObjF3(p)-ObjF3)<0 )  | ...
                                ((ObjF1(p)-ObjF1)<0  & (ObjF2(p)-ObjF2)==0 & (ObjF3(p)-ObjF3)<0 )  | ...
                                ((ObjF1(p)-ObjF1)<0  & (ObjF2(p)-ObjF2)==0 & (ObjF3(p)-ObjF3)==0 ) | ...
                                ((ObjF1(p)-ObjF1)<0  & (ObjF2(p)-ObjF2)<0  & (ObjF3(p)-ObjF3)==0 ) | ...
                                ((ObjF1(p)-ObjF1)==0 & (ObjF2(p)-ObjF2)<0  & (ObjF3(p)-ObjF3)==0 ));

            
            n(p) = length(find(((ObjF1(p)-ObjF1)>0  & (ObjF2(p)-ObjF2)>0  & (ObjF3(p)-ObjF3)>0 )  | ...
                               ((ObjF1(p)-ObjF1)==0 & (ObjF2(p)-ObjF2)>0  & (ObjF3(p)-ObjF3)>0 )  | ...
                               ((ObjF1(p)-ObjF1)==0 & (ObjF2(p)-ObjF2)==0 & (ObjF3(p)-ObjF3)>0 )  | ...
                               ((ObjF1(p)-ObjF1)>0  & (ObjF2(p)-ObjF2)==0 & (ObjF3(p)-ObjF3)>0 )  | ...
                               ((ObjF1(p)-ObjF1)>0  & (ObjF2(p)-ObjF2)==0 & (ObjF3(p)-ObjF3)==0 ) | ...
                               ((ObjF1(p)-ObjF1)>0  & (ObjF2(p)-ObjF2)>0  & (ObjF3(p)-ObjF3)==0 ) | ...
                               ((ObjF1(p)-ObjF1)==0 & (ObjF2(p)-ObjF2)>0  & (ObjF3(p)-ObjF3)==0 )));
        end
        
        % First front
        front(1).fr = find(n==0);

        % Subsequent fronts
        while (~isempty(front(frontCount).fr))      
            lastFront = front(frontCount).fr;    
            n(lastFront) = inf;                  
            feasSolutions(lastFront,G+F+1) = frontCount; 
            frontCount = frontCount + 1;            
            front(frontCount).fr=[];            
           for i = 1:length(lastFront)        
                tempf = dominationSet(lastFront(i)).sp; 
                n(tempf)=n(tempf)-1;            
           end

                q = find(n==0);                                  
                front(frontCount).fr = [front(frontCount).fr q];  

        end 
          
        % Riordino i cromosomi in ordine secondo il fronte a cui appartengo
        chromosome_sorted = sortrows(feasSolutions,G+F+1); 

        %% Calcolo della Crowding Distance di ogni individuo
        rowsindex = 1;

        for i = 1:length(front) - 1         % Per i che va da 1 al numero di fronti - 1 (non considero l'ultimo fronte che e' vuoto) 
            length_f = length(front(i).fr); % Ricavo la lunghezza del fronte considerato

             if length_f > 2                % Se la lunghezza del fronte considerato e' > 2

              sorted_index_f1=[];           % Iniziallizzo le variabili
              sorted_index_f2=[];
              sorted_f1=[];
              sorted_f2=[];

            % Ordino gli individui del fronte considerato in base al loro valore di ciascuna
            % funzione obiettivo, prima per f1 e poi per f2.
            % sorted_f1 contiene i valori della funzione obiettivo f1 assunti dai
            % vari individui del fronte considerato ordinati dal minore al maggiore
            % (ugualmente in sorted_f2 per la funzione obiettivo f2).
            % sorted_index_f1 contiene gli indici degli individui del fronte
            % considerato ordinati da quello che ha valore minore di f1 a quello che ha
            % valore maggiore. (analogamente in sorted_index_f2 per la funzione f2).
            [sorted_f1, sorted_index_f1] = sortrows(chromosome_sorted(rowsindex:(rowsindex + length_f - 1),G+1));
            [sorted_f2, sorted_index_f2] = sortrows(chromosome_sorted(rowsindex:(rowsindex + length_f - 1),G+2));
            [sorted_f3, sorted_index_f3] = sortrows(chromosome_sorted(rowsindex:(rowsindex + length_f - 1),G+3));

            % Ricavo il valore minimo e massimo della funzione f1 assunto dagli
            % individui del fronte considerato
            f1min = chromosome_sorted(sorted_index_f1(1) + rowsindex - 1,G+1);
            f1max = chromosome_sorted(sorted_index_f1(end) + rowsindex - 1,G+1);

            % Creo una nuova colonna nella matrice chromosome_sorted.
            % Tale colonna indica la CD relativa a f1 di ciascun individuo,
            % pongo il valore infinito in corrispondenza degli individui che hanno
            % valore minimo e massimo della funzione obiettivo f1
            chromosome_sorted(sorted_index_f1(1)+rowsindex-1,G+F+2) = inf;
            chromosome_sorted(sorted_index_f1(end)+rowsindex-1,G+F+2) = inf;

            % Faccio l'analogo per f2:
            % Ricavo il valore minimo e massimo della funzione f2 assunto dagli
            % individui del fronte considerato
            f2min = chromosome_sorted(sorted_index_f2(1) + rowsindex - 1,G+2);
            f2max = chromosome_sorted(sorted_index_f2(end) + rowsindex - 1,G+2);

            % Creo una nuova colonna nella matrice chromosome_sorted.
            % Tale colonna indica la CD relativa a f2 di ciascun individuo,
            % pongo il valore infinito in corrispondenza degli individui che hanno
            % valore minimo e massimo della funzione obiettivo f2
            chromosome_sorted(sorted_index_f2(1) + rowsindex - 1,G+F+3) = inf;
            chromosome_sorted(sorted_index_f2(end) + rowsindex - 1,G+F+3) = inf;
            
            % Faccio l'analogo per f3:
            % Ricavo il valore minimo e massimo della funzione f3 assunto dagli
            % individui del fronte considerato
            f3min = chromosome_sorted(sorted_index_f3(1) + rowsindex - 1,G+3);
            f3max = chromosome_sorted(sorted_index_f3(end) + rowsindex - 1,G+3);

            % Creo una nuova colonna nella matrice chromosome_sorted.
            % Tale colonna indica la CD relativa a f2 di ciascun individuo,
            % pongo il valore infinito in corrispondenza degli individui che hanno
            % valore minimo e massimo della funzione obiettivo f2
            chromosome_sorted(sorted_index_f3(1) + rowsindex - 1,G+F+4) = inf;
            chromosome_sorted(sorted_index_f3(end) + rowsindex - 1,G+F+4) = inf;

            % Calcolo la CD degli altri elementi del fronte
             for j = 2:length(front(i).fr) - 1
                 % Se il valore minimo e massimo di una delle 3 funz ob coincidono
                 % Pongo infinito i valori sulle 2 CD per ogni
                 % individuo del fronte.
                 if  (f1max - f1min == 0) || (f2max - f2min == 0) || (f3max - f3min == 0)
                     chromosome_sorted(sorted_index_f1(j) + rowsindex - 1,G+F+2) = inf;
                     chromosome_sorted(sorted_index_f2(j) + rowsindex - 1,G+F+3) = inf;
                     chromosome_sorted(sorted_index_f3(j) + rowsindex - 1,G+F+4) = inf;
                 else                                               % Altrimenti
     % Calcolo il valore della crowding distance secondo la formula usuale, per ciascuna delle due funzioni:
     chromosome_sorted(sorted_index_f1(j) + rowsindex - 1,G+F+2) = ...
     (chromosome_sorted(sorted_index_f1(j+1) + rowsindex - 1,G+1) - ...
     chromosome_sorted(sorted_index_f1(j-1) + rowsindex - 1,G+1))/(f1max - f1min);
     chromosome_sorted(sorted_index_f2(j) + rowsindex - 1,G+F+3) = ...
     (chromosome_sorted(sorted_index_f2(j+1) + rowsindex - 1,G+2) - ...
     chromosome_sorted(sorted_index_f2(j-1) + rowsindex - 1,G+2))/(f2max - f2min);
     chromosome_sorted(sorted_index_f3(j) + rowsindex - 1,G+F+4) = ...
     (chromosome_sorted(sorted_index_f3(j+1) + rowsindex - 1,G+3) - ...
     chromosome_sorted(sorted_index_f3(j-1) + rowsindex - 1,G+3))/(f3max - f3min);
                 end
             end

             else % Se nel fronte ci sono <= 2 individui pongo infinito i loro valori delle CD
                chromosome_sorted(rowsindex:(rowsindex + length_f - 1),G+F+2:G+F+4) = inf;
              end
         rowsindex = rowsindex + length_f; % Aumento l'indice di riga per i successivi fronti
        end
        % Aggiungo una colonna in cui sommo i valori della crowding distance
        % relativi alle due funz obiettivo, ottenendo la CD finale
        % dell'individuo.
        chromosome_sorted(:,G+F+5) = sum(chromosome_sorted(:,G+F+2:G+F+4),2); 

    % L'output finale e' la matrice dei cromosomi in cui gli individui sono
    % ordinati in base al fronte in cui stanno, e nell'ultima colonna e'
    % esplicitata la crowding distance di ciascun individuo.
    populationSorted1 = [chromosome_sorted(:,1:G+F) zeros(feasSize,1) chromosome_sorted(:,G+F+1) chromosome_sorted(:,G+F+5)];
    end
    % Se ci sono soluzioni non ammissibili assegno ad esse fronte crescente
    % al crescere dell'errore, a partire dall'ultimo fronte delle soluzioni
    % ammissibili. Assegno inoltre CD Inf.
    if problemType=='u' || problemType=='m'       
        infeasiblePop = sortrows(infSolutions,G+F+1);
        infeasiblePop = [infeasiblePop(:,1:G+F+1) (frontCount:frontCount-1+size(infeasiblePop,1))' inf*(ones(size(infeasiblePop,1),1))];
        for kk = (size(front,2)):(size(front,2))+(length(infSolutions))-1
         front(kk).fr= feasSize+1;
        end
    end
populationSorted = [populationSorted1;infeasiblePop];
