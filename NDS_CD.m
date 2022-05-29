function [populationSorted] = NDS_CD(population)

    global G F problemType

    %% Iniziallizzazione delle variabili
    dominationSet.sp = [];
    frontCount = 1;
    populationSorted1 = [];
    infeasiblePop = [];
    front.fr=[];
    %% Feasible solutions identification
    if all(population(:,G+F+1)==0)              
        problemType = 'f';                   
        feasSolutions = population(:,1:G+F);      
        feasSize = size(feasSolutions,1);      
    elseif all(population(:,G+F+1)~=0)       
        problemType = 'u';                    
        feasSize = 0;                     
        infSolutions = population;    
    else
        problemType = 'm';
        feasIndex = find(population(:,G+F+1)==0);   
        feasSolutions = population(feasIndex,1:G+F);   
        feasSize = size(feasSolutions,1);                
        infeasIndex = find(population(:,G+F+1)~=0);   
        infSolutions = population(infeasIndex,1:G+F+1); 
    end
    
    
    %%
    if problemType == 'f' || problemType == 'm'   
        
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
       
        front(1).fr = find(n==0);

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
           
        sortedPop = sortrows(feasSolutions,G+F+1); 

        %% CD
        rowsIndex = 1;

        for i = 1:length(front) - 1   
            lengthFront = length(front(i).fr); 

             if lengthFront > 2             
              
            [~, sortedIndexObjF1] = sortrows(sortedPop(rowsIndex:(rowsIndex + lengthFront - 1),G+1));
            [~, sortedIndexObjF2] = sortrows(sortedPop(rowsIndex:(rowsIndex + lengthFront - 1),G+2));
            [~, sortedIndexObjF3] = sortrows(sortedPop(rowsIndex:(rowsIndex + lengthFront - 1),G+3));

            ObjF1min = sortedPop(sortedIndexObjF1(1) + rowsIndex - 1,G+1);
            ObjF1max = sortedPop(sortedIndexObjF1(end) + rowsIndex - 1,G+1);
            
            sortedPop(sortedIndexObjF1(1)+rowsIndex-1,G+F+2) = inf;
            sortedPop(sortedIndexObjF1(end)+rowsIndex-1,G+F+2) = inf;

            
            ObjF2min = sortedPop(sortedIndexObjF2(1) + rowsIndex - 1,G+2);
            ObjF2max = sortedPop(sortedIndexObjF2(end) + rowsIndex - 1,G+2);
            
            sortedPop(sortedIndexObjF2(1) + rowsIndex - 1,G+F+3) = inf;
            sortedPop(sortedIndexObjF2(end) + rowsIndex - 1,G+F+3) = inf;
            
            
            ObjF3min = sortedPop(sortedIndexObjF3(1) + rowsIndex - 1,G+3);
            ObjF3max = sortedPop(sortedIndexObjF3(end) + rowsIndex - 1,G+3);

            sortedPop(sortedIndexObjF3(1) + rowsIndex - 1,G+F+4) = inf;
            sortedPop(sortedIndexObjF3(end) + rowsIndex - 1,G+F+4) = inf;

            
             for j = 2:length(front(i).fr) - 1
                 if  (ObjF1max - ObjF1min == 0) || (ObjF2max - ObjF2min == 0) || (ObjF3max - ObjF3min == 0)
                     sortedPop(sortedIndexObjF1(j) + rowsIndex - 1,G+F+2) = inf;
                     sortedPop(sortedIndexObjF2(j) + rowsIndex - 1,G+F+3) = inf;
                     sortedPop(sortedIndexObjF3(j) + rowsIndex - 1,G+F+4) = inf;
                 else                                              
     sortedPop(sortedIndexObjF1(j) + rowsIndex - 1,G+F+2) = ...
     (sortedPop(sortedIndexObjF1(j+1) + rowsIndex - 1,G+1) - ...
     sortedPop(sortedIndexObjF1(j-1) + rowsIndex - 1,G+1))/(ObjF1max - ObjF1min);
     sortedPop(sortedIndexObjF2(j) + rowsIndex - 1,G+F+3) = ...
     (sortedPop(sortedIndexObjF2(j+1) + rowsIndex - 1,G+2) - ...
     sortedPop(sortedIndexObjF2(j-1) + rowsIndex - 1,G+2))/(ObjF2max - ObjF2min);
     sortedPop(sortedIndexObjF3(j) + rowsIndex - 1,G+F+4) = ...
     (sortedPop(sortedIndexObjF3(j+1) + rowsIndex - 1,G+3) - ...
     sortedPop(sortedIndexObjF3(j-1) + rowsIndex - 1,G+3))/(ObjF3max - ObjF3min);
                 end
             end

             else 
                sortedPop(rowsIndex:(rowsIndex + lengthFront - 1),G+F+2:G+F+4) = inf;
              end
         rowsIndex = rowsIndex + lengthFront; 
        end
        sortedPop(:,G+F+5) = sum(sortedPop(:,G+F+2:G+F+4),2); 

    
    populationSorted1 = [sortedPop(:,1:G+F) zeros(feasSize,1) sortedPop(:,G+F+1) sortedPop(:,G+F+5)];
    end
    if problemType=='u' || problemType=='m'       
        infeasiblePop = sortrows(infSolutions,G+F+1);
        infeasiblePop = [infeasiblePop(:,1:G+F+1) (frontCount:frontCount-1+size(infeasiblePop,1))' inf*(ones(size(infeasiblePop,1),1))];
        for kk = (size(front,2)):(size(front,2))+(length(infSolutions))-1
         front(kk).fr= feasSize+1;
        end
    end
populationSorted = [populationSorted1;infeasiblePop];