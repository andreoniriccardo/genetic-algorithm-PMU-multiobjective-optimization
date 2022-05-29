clear variables
clear global
close all
clc
%% Initialization of variables
% NSGA II parameters
global G F E N mutProb
% G -> number of genes (optimization variables);
% F -> number of objective functions;
% E -> number of constraints;
% N -> population size (number of individuals);
% finalGen -> number of generations;
% mutProb -> mutation probability;

% Definition of the power system
test_case = case141;

% Power system parameters
global Incidence Gbus Bbus V_R V_I Iinj_R Iinj_I x_nom n_meas gen max_f2
global R_V_R R_V_I R_Iinj_R R_Iinj_I W_V_R W_V_I W_Iinj_R W_Iinj_I std_V_R std_V_I std_Iinj_R std_Iinj_I 
global ZIbus Gbus_ZI Bbus_ZI W_ZI_IR W_ZI_II R_ZI_IR R_ZI_II omega nbranch Branch_connectivity
global C_line g_line b_line Ibranch_R Ibranch_I R_Ibranch_R R_Ibranch_I
global W_Ibranch_R W_Ibranch_I std_Ibranch_R std_Ibranch_I cost n_sim_sens delta_H
global Gbus_noZI Bbus_noZI nnzGBbus_noZI nnzGBbus_ZI nnzgbline H_max_nominal R_max
global nnzBbus_noZI
global M11_f3 M12_f3 M21_f3 M22_f3 Gbus_NoZI_f3 Bbus_NoZI_f3 Gbus_ZI_f3 Bbus_ZI_f3
[Incidence, G, x_nom, Gbus, Bbus, V_R, V_I, Iinj_R, Iinj_I, ...
 R_V_R, R_V_I, R_Iinj_R, R_Iinj_I, W_V_R, W_V_I, W_Iinj_R, W_Iinj_I, ...
 std_V_R, std_V_I, std_Iinj_R, std_Iinj_I, ZIbus, Gbus_ZI, Bbus_ZI, ...
 W_ZI_IR, W_ZI_II, R_ZI_IR, R_ZI_II, omega, nbranch, Branch_connectivity, ...
 C_line, g_line, b_line, Ibranch_R, Ibranch_I, R_Ibranch_R, R_Ibranch_I, ...
 W_Ibranch_R, W_Ibranch_I, std_Ibranch_R, std_Ibranch_I, cost] = pmu_data(test_case);
save('cost_141new.mat','cost')
% Setting parameters values
n_meas = 100; % Number of measures provided by each PMU
F = 3;
% Ci sono gli G vincoli tradizionali più nbranch*G vincoli sulle contingencies di linea 
% più G*G vincoli sulle contingencies delle pmu
E = G+G*nbranch+G*G;    
N = 1000;
mutProb = 0.1; % Probabilita' di mutazione
finalGen = 120;
gen = 1;
% Evaluation of constraints compliance and objective functions
obj_fun_eval = 'pmu_evalZI_correnteLinee';

% Memory matrix initiallization
% - Genome                      (columns 1:G)
% - Objective function values   (columns G+1:G+F)
% - Normalized error            (column  G+F+1)
% - Front                       (column  G+F+2)
% - CD                          (column  G+F+3)
% - Generation                  (column  G+F+4)
Memory = zeros(N*finalGen, G+F+1+2+1);
max_f2 = 1e2*ones(finalGen,1);
max_f2(1:finalGen/4) = 1e3;
%% Creazione della matrice delta_H per la misura della sensibilità
% Calcolo la H massima (nel caso di piazzamento di tutte le PMU)

%{
Gbus_noZI = Gbus;
Gbus_noZI(ZIbus,:) = [];
Bbus_noZI = Bbus;
Bbus_noZI(ZIbus,:) = [];

H_max_nominal = [diag(ones(G,1)), zeros(G,G);
                 zeros(G,G),      diag(ones(G,1));
                 Gbus_noZI,      -Bbus_noZI;
                 Bbus_noZI,       Gbus_noZI;
                 Gbus_ZI,        -Bbus_ZI;
                 Bbus_ZI,         Gbus_ZI;
                 g_line,         -b_line;
                 b_line,          g_line];                
% Calcolo la corrispondente matrice R
R_Iinj_R_noZI = R_Iinj_R;
R_Iinj_R_noZI(ZIbus,:) = [];
R_Iinj_I_noZI = R_Iinj_I;
R_Iinj_I_noZI(ZIbus,:) = [];

diagR = [R_V_R; R_V_I; R_Iinj_I_noZI; R_Iinj_I_noZI; R_ZI_IR; R_ZI_II; R_Ibranch_R; R_Ibranch_I];
R_max = diag(diagR);

nnzGBbus_noZI = nnz(Gbus_noZI);
nnzBbus_noZI = nnz(Bbus_noZI);
nnzGBbus_ZI = nnz(Gbus_ZI);
nnzgbline = nnz(g_line);
nnzbline = nnz(b_line);

sizex = nnzGBbus_noZI + nnzGBbus_ZI + nnzgbline;
%x0 = ones(sizex,1);
%x = fminsearch(@funmaxS,x0)
lb = 0.9*ones(sizex,1);
ub = 1.1*ones(sizex,1);
x_ = particleswarm(@funmaxS,sizex,lb,ub) 

M11 = diag(ones(G,1));
M12 = zeros(G,G);
M21 = zeros(G,G);
M22 = diag(ones(G,1));

xGBbus_noZI = x(1:nnzGBbus_noZI);
xGBbus_ZI = x(nnzGBbus_noZI+1:nnzGBbus_noZI+nnzGBbus_ZI);
xgbline = x(nnzGBbus_noZI+nnzGBbus_ZI+1:nnzGBbus_noZI+nnzGBbus_ZI+nnzgbline);

[row_Gbus_noZI, col_Gbus_noZI] = find(Gbus_noZI);
MGbus_noZI = accumarray([row_Gbus_noZI(:),col_Gbus_noZI(:)],xGBbus_noZI(:));
MBbus_noZI = accumarray([row_Gbus_noZI(:),col_Gbus_noZI(:)],xGBbus_noZI(:));
% Correggo elementi sbagliati di MB solo per rete a 141 nodi
MBbus_noZI(50,86)=1;
MBbus_noZI(49,87)=1;
MBbus_noZI(50,87)=1;

[row_Gbus_ZI, col_Gbus_ZI] = find(Gbus_ZI);
MGbus_ZI = accumarray([row_Gbus_ZI(:),col_Gbus_ZI(:)],xGBbus_ZI(:));
MBbus_ZI = accumarray([row_Gbus_ZI(:),col_Gbus_ZI(:)],xGBbus_ZI(:));
[abc, ncol_MGbus_ZI] = size(MGbus_ZI);
[abc, ncol_MBbus_ZI] = size(MBbus_ZI);
if ncol_MGbus_ZI < G
    diff = G-ncol_MGbus_ZI;
    MGbus_ZI(:,ncol_MGbus_ZI+1:G)=zeros(numel(ZIbus,diff));
end
if ncol_MBbus_ZI < G
    diff = G-ncol_MBbus_ZI;
    MBbus_ZI(:,ncol_MBbus_ZI+1:G)=zeros(numel(ZIbus,diff));
end

[row_gline, col_gline] = find(g_line);
Mgline = accumarray([row_gline(:),col_gline(:)],xgbline(:));
Mbline = accumarray([row_gline(:),col_gline(:)],xgbline(:));
% Correggo elementi sbagliati di Mb solo per rete a 141 nodi
Mbline(86,86) = 1;
Mbline(86,87) = 1;


Tol = [M11              M12;
       M21              M22;
       MGbus_noZI      -MBbus_noZI;
       MBbus_noZI       MGbus_noZI;
       MGbus_ZI        -MBbus_ZI;
       MBbus_ZI         MGbus_ZI;
       Mgline          -Mbline
       Mbline           Mgline];

delta_H = H_max_nominal.*Tol;
save('deltaH_141nuovo_def.mat','delta_H')
%}
%{
n_sim_sens = 100; % Numero di simulazioni 
% Calcolo la H massima (nel caso di piazzamento di tutte le PMU)

Gbus_noZI = Gbus;
Gbus_noZI(ZIbus,:) = [];
Bbus_noZI = Bbus;
Bbus_noZI(ZIbus,:) = [];

H_max_nominal = [diag(ones(G,1)), zeros(G,G);
                 zeros(G,G),      diag(ones(G,1));
                 Gbus_noZI,      -Bbus_noZI;
                 Bbus_noZI,       Gbus_noZI;
                 Gbus_ZI,        -Bbus_ZI;
                 Bbus_ZI,         Gbus_ZI;
                 g_line,         -b_line;
                 b_line,          g_line];
                 
% Iniziallizzo la matrice delta_H che raccoglie n_sim_sens matrici H max
% con i valori randomizzati
%delta_H = zeros(n_sim_sens*(G+G+G+G+numel(ZIbus)+numel(ZIbus) + 2*nbranch),2*G);
for i = 1:n_sim_sens
    H_rand = (ones(size(H_max_nominal))+0.2/3*randn(size(H_max_nominal)));
    H_rand(1:G,1:G) = diag(ones(G,1));
    H_rand(1:G,G+1:2*G) = zeros(G,G);
    H_rand(G+1:2*G,1:G) = zeros(G,G);
    H_rand(G+1:2*G,G+1:2*G) = diag(ones(G,1));
    delta_H(1:G+G+2*G+ 2*nbranch,1:2*G,i)= H_max_nominal.*H_rand;
end
%}
C = load('deltaH_141nuovo_def.mat');
delta_H = C.delta_H;

Gbus_NoZI_f3 = delta_H(2*G+1:2*G+(G-numel(ZIbus)),1:G);
Bbus_NoZI_f3 = delta_H(2*G+(G-numel(ZIbus))+1:2*G+2*(G-numel(ZIbus)),1:G);

Gbus_ZI_f3 = delta_H(2*G+2*(G-numel(ZIbus))+1:2*G+2*(G-numel(ZIbus))+numel(ZIbus),1:G);
Bbus_ZI_f3 = delta_H(2*G+2*(G-numel(ZIbus))+numel(ZIbus)+1:2*G+2*(G-numel(ZIbus))+2*numel(ZIbus),1:G);

M11_f3 = diag(ones(G,1));
M12_f3 = zeros(G,G);
M21_f3 = zeros(G,G);
M22_f3 = diag(ones(G,1));


%% Generation of the initial population
%{
xLow = zeros(1,G);       % Lower limit of genome
xUpp = ones(1,G);        % Upper limit of genome

popLow = repmat(xLow, N,1);    % Lower population limit
popUpp = repmat(xUpp, N,1);    % Upper population limit

% Random generation of the initial populaiton
popInit = popLow + ((popUpp - popLow).*rand(N,G));
popInit = round(popInit);
%}
popInit = initialPopOptStimatoreN();

% Evaluation of objective functions and constraint compliance of initial
% population
objFuns = zeros(N,F);
errors  = zeros(N,E);
for i = 1:N
   [objFuns(i,:), errors(i,:)] = feval(obj_fun_eval, popInit(i,:)); 
   disp({'Evaluating objective functions of the initial population',i});
end
% Error normalization
errorNorm = normalize(errors);

% Initial population
initialPopulation = [popInit objFuns errorNorm];

% Initial population non-dominated sorting and CD computation
[population] = NDS_CD(initialPopulation);

h = figure(1);
%% Generations

for gen = 1:finalGen

    % Parents selection
    selectedParents = tournamentSelection(population);

    % Generazione della popolazione figlia Q_t trammite crossover
    offspring  = geneticOperator(selectedParents(:,1:G));

    % Offspring population evaluation
    offObjFuns = zeros(N,F);
    offErrors  = zeros(N,E);
    for ii = 1:N
        [offObjFuns(ii,:), offErrors(ii,:)] = feval(obj_fun_eval, offspring(ii,:));
        disp({'Evaluating objective functions of the generation n. ',gen,ii});
    end
    % Normalized errors
    errorNorm = normalize(offErrors);
    
    % Offspring populaiton
    offspringPopulation = [offspring offObjFuns errorNorm];

    % Itermediate population
    interPopulation = [population(:,1:G+F+1); offspringPopulation(:,1:G+F+1)];

    % Intermediate population sorting
    [interPopSorted] = NDS_CD(interPopulation);

    % Selection of N individuals
    newPopulation = replace(interPopSorted);
    population = newPopulation;
    
    %% Plot
    hold off;        
    figure(1)
    plot(newPopulation(:,G+1),newPopulation(:,G+2),'o')
    title({'PMU placement optimization'; ['Generation ', num2str(gen)]})
    xlabel('Number of Channels')
    ylabel('Measurement Uncertainty')
    grid on
    drawnow
    figure(2)
    plot(newPopulation(:,G+1),newPopulation(:,G+3),'o')
    title({'PMU placement optimization'; ['Generation ', num2str(gen)]})
    xlabel('Number of Channels')
    ylabel('Sensitivity')
    grid on
    drawnow
    %% Storing the previous generations data
    Memory((gen-1)*N+1:gen*N,1:G+F+1+2) = newPopulation;
    Memory((gen-1)*N+1:gen*N,G+F+1+2+1)=gen*ones(N,1);   
    
end
save('141bus_1000_010_0100_120_01.mat','newPopulation')
save('Memory_141bus_1000_010_0100_120_01.mat','Memory')