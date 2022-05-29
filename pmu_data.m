function [Inc, nbus, x_nom, Gbus, Bbus, V_R, V_I, Iinj_R, Iinj_I, ...
          R_V_R, R_V_I, R_Iinj_R, R_Iinj_I, W_V_R, W_V_I, W_Iinj_R, W_Iinj_I, ...
          std_V_R, std_V_I, std_Iinj_R, std_Iinj_I, ZIbus, Gbus_ZI, Bbus_ZI, ...
          W_ZI_IR, W_ZI_II, R_ZI_IR, R_ZI_II, omega, nbranch, Branch_connectivity, ...
          C_line, g_line, b_line, Ibranch_R, Ibranch_I, R_Ibranch_R, R_Ibranch_I, ...
          W_Ibranch_R, W_Ibranch_I, std_Ibranch_R, std_Ibranch_I, cost] = pmu_data(x)
%% Caricamento della rete ed estrapolazione info preliminari
% Caricamento del test case specificato in input
mpc = x;
% Definizione del numero di bus e di linee della rete
nbus_ = size(mpc.bus);
nbus  = nbus_(1);
nbranch_ = size(mpc.branch);
nbranch  = nbranch_(1);
% Ricavo i bus ZI
ZIcount = 1;
for i = 2:nbus
    if mpc.bus(i,3) == 0 && mpc.bus(i,4) == 0 && mpc.bus(i,2) ~= 2        
       ZIbus(ZIcount) = i; 
       ZIcount = ZIcount + 1;
    end
end
if ZIcount == 1
    ZIbus = [];
end
%% Estrapolazione matrici topologiche della rete
% Creazione della matrice di incidenza della linea
Ainc = makeIncidence(mpc);
% La matrice bus indica con quali linee è comunicante ciascun bus.
% L'elemento bus(i, j) se diverso da zero indica che il bus j è
% collegato alla linea i.
bus = zeros(nbus,nbranch);
for i = 1:nbus
    ind = 1;
    for j = 1:nbranch        
        if Ainc(j,i) == 1
            bus(i,ind) = j;
            ind = ind + 1;
        elseif Ainc(j,i) == -1
            bus(i,ind) = j;
            ind = ind + 1;
        else
            bus(i,ind) = 0;
            ind = ind + 1;
        end
    end
end

% Creazione della matrice di incidenza Inc
Inc = diag(ones(nbus,1));
% Trovo i bus collegati alla medesima linea ii, tali bus sono adiacenti.
for ii = 1:nbranch           
    temp = find(bus(:,ii));
    Inc(temp(1),temp(2)) = 1;
    Inc(temp(2),temp(1)) = 1;
end

% Creazione della matrice di connettività delle linee:
% la riga i mi indica quali bus sono connessi dalla linea i
Branch_connectivity = zeros(nbranch,2);
for i = 1:nbranch
    ind = 1;
    for j = 1:nbus
        if Ainc(i,j) == 1 || Ainc(i,j) == -1
            Branch_connectivity(i,ind) = j;
            ind = ind + 1;
        end
    end
end

% Creo matrice della connettività delle linee C_line, simile alla precedente ma
% scritta in modo diverso:
% Nella riga i-esima della matrice C_line tutti gli elementi sono 0 tranne
% quelli occupanti le colonne corrispondenti ai bus collegati dalla linea
% i-esima. Quindi la matrice C_line ha dimensione n_branh X n_bus
C_line_provv = makeIncidence(mpc);
C_line = zeros(nbranch,nbus);
% Rappresentazione in forma matriciale standard
for i = 1:nbranch
   for j = 1:nbus
      C_line(i,j) = C_line_provv(i,j);
   end    
end

% Ricavo vettore dei costi, indicante quanti PMU vanno piazzati su un nodo.
% Nel nodo i vanno piazzati cost_i = (N_i + 2)/2 PMU arrotondato per
% eccesso.
N_inc_branch = zeros(1,nbus); % Numero di linee incidenti ad ogni nodo
for i = 1:nbus
    N_inc_branch(i) = nnz(Inc(i,:)) - 1; % Meno 1 perchè non considero la diagonale della matrice delle incidenze    
end
cost = zeros(1,nbus);
for i = 1:nbus
    %cost(i) = ceil((N_inc_branch(i) + 2)/2); % ceil per arrotondare per eccesso    
    cost(i) = N_inc_branch(i) + 2;
end




% Ricavo matrice delle ammettenze Ybus
% e le matrici che collegano corrente sulle linee con le tensioni ai nodi
% Yf e Yt (correndi DA e VERSO la fine delle linee)
[Ybus_, Yf_, Yt_] = makeYbus(mpc);
Ybus = zeros(nbus,nbus);
% Rappresentazione in forma matriciale standard
for i = 1:nbus
   for j = 1:nbus
      Ybus(i,j) = Ybus_(i,j);
   end    
end
% Suddivisione della matrice Ybus in parte reale ed immaginaria
Gbus = real(Ybus); % Parte reale di Ybus
Bbus = imag(Ybus); % Parte immaginaria di Ybus

Yf = zeros(nbranch,nbus);
for i = 1:nbranch
   for j = 1:nbus
      Yf(i,j) = Yf_(i,j);
   end    
end
Yt = zeros(nbranch,nbus);
for i = 1:nbranch
   for j = 1:nbus
      Yt(i,j) = Yt_(i,j);
   end    
end
% Suddivisione della matrice Yf in parte reale ed immaginaria
g_line = real(Yf);
b_line = imag(Yf);

%% Estrapolazione dati relativi a ZI
% Individuo nelle matrici Gbus e Bbus le righe corrispondenti agli ZI
% bus
if ZIcount > 1 % Se c'è qualche ZI bus
    Gbus_ZI = zeros(numel(ZIbus), nbus);
    Bbus_ZI = zeros(numel(ZIbus), nbus);

    for i = 1:numel(ZIbus)
       Gbus_ZI(i,:) = Gbus(ZIbus(i),:);
       Bbus_ZI(i,:) = Bbus(ZIbus(i),:);        
    end

    % Cardinalita' regioni ZI
    omega = zeros(numel(ZIbus),1);
    for i = 1:numel(ZIbus)
       omega(i) = nnz(Inc(ZIbus(i),:)); 
    end
else
    Gbus_ZI = [];
    Bbus_ZI = [];
    omega = [];
end
%% Risoluzione della rete ed estrapolazione del ground truth
% Risoluzione della rete
results = runpf(mpc);

% Magnitudine della tensione ai bus (in p.u.),
% indicate nella colonna 8 di result.bus:
V_mag = results.bus(:,8);
% Angolo della tensione  ai bus (in deg),
% indicati nella colonna 9 di result.bus:
V_ang = results.bus(:,9);
V_ang_rad = deg2rad(V_ang);
    
% Grandezze in coordinate rettangolari
% Tensioni nominali in coordinate rettangolari (componente reale ed
% immaginaria)
V_R = V_mag.*cos(V_ang_rad);
V_I = V_mag.*sin(V_ang_rad);
% Stato nominale
x_nom = [V_R; V_I];

% Correnti iniettate nette ai bus in coordinate rettangolari
Iinj_R = zeros(nbus, 1);
Iinj_I = zeros(nbus, 1);
for i = 1:nbus    
    Iinj_R(i) = dot(V_R,Gbus(i,:)) - dot(V_I,Bbus(i,:));
    Iinj_I(i) = dot(V_R,Bbus(i,:)) + dot(V_I,Gbus(i,:));
end    
% Correnti iniettate nette ai bus in coordinate polari
Iinj_mag = zeros(nbus, 1);
Iinj_ang = zeros(nbus, 1);
for i = 1:nbus
    Iinj_mag(i) = sqrt(Iinj_R(i)^2 + Iinj_I(i)^2);
    Iinj_ang(i) = atan2(Iinj_I(i), Iinj_R(i));    
end

% Correnti di linea in coordinate rettangolari
Ibranch_R = zeros(nbranch,1);
Ibranch_I = zeros(nbranch,1);
for i = 1:nbranch
    Ibranch_R(i) = dot(g_line(i,:), V_R) - dot(b_line(i,:), V_I);
    Ibranch_I(i) = dot(b_line(i,:), V_R) + dot(g_line(i,:), V_I);    
end
% Correnti di linea in coordinate polari
Ibranch_mag = zeros(nbranch, 1);
Ibranch_ang = zeros(nbranch, 1);
for i = 1:nbranch
    Ibranch_mag(i) = sqrt(Ibranch_R(i)^2 + Ibranch_I(i)^2);
    Ibranch_ang(i) = atan2(Ibranch_I(i), Ibranch_R(i));
end

%{
for i = 1:numel(ZIbus)
    Iinj_R(ZIbus(i)) = NaN;
    Iinj_I(ZIbus(i)) = NaN;
end
%}
%% Varianze
% Deviazioni standard coordinate polari
sigma_V = 0.0033*abs(V_mag);
sigma_thetaV = 0.0033*ones(nbus,1);
sigma_Iinj = 0.0033*abs(Iinj_mag);
sigma_thetaIinj = 0.0033*ones(nbus,1);
sigma_Ibranch =0.0033*abs(Ibranch_mag);
sigma_thetaIbranch = 0.0033*ones(nbranch,1);

% Conversione delle varianze delle grandezze da coordinate polari a
% coordinate rettangolari
%{
for i = 1:nbus
   sigma2_V_R(i) = (V_mag(i)^2)*((sin(V_ang_rad(i)))^2)*(sigma_thetaV(i)^2) +...
                    (sigma_V(i)^2)*((cos(V_ang_rad(i)))^2) - ...
                    V_mag(i)*sigma_V(i)*sigma_thetaV(i)*sin(2*V_ang_rad(i));
   sigma2_V_I(i) = (V_mag(i)^2)*((cos(V_ang_rad(i)))^2)*(sigma_thetaV(i)^2) +...
                    (sigma_V(i)^2)*((sin(V_ang_rad(i)))^2) + ...
                    V_mag(i)*sigma_V(i)*sigma_thetaV(i)*sin(2*V_ang_rad(i));
                
  sigma2_Iinj_R(i) = (Iinj_mag(i)^2)*((sin(Iinj_ang(i)))^2)*(sigma_thetaIinj(i)^2) +...
                    (sigma_Iinj(i)^2)*((cos(Iinj_ang(i)))^2) - ...
                    Iinj_mag(i)*sigma_Iinj(i)*sigma_thetaIinj(i)*sin(2*Iinj_ang(i));
  sigma2_Iinj_I(i) = (Iinj_mag(i)^2)*((cos(Iinj_ang(i)))^2)*(sigma_thetaIinj(i)^2) +...
                    (sigma_Iinj(i)^2)*((sin(Iinj_ang(i)))^2) + ...
                    Iinj_mag(i)*sigma_Iinj(i)*sigma_thetaIinj(i)*sin(2*Iinj_ang(i));
end
for i = 1:nbranch
  sigma2_Ibranch_R(i) = (Ibranch_mag(i)^2)*((sin(Ibranch_ang(i)))^2)*(sigma_thetaIbranch(i)^2) +...
                    (sigma_Ibranch(i)^2)*((cos(Ibranch_ang(i)))^2) - ...
                    Ibranch_mag(i)*sigma_Ibranch(i)*sigma_thetaIbranch(i)*sin(2*Ibranch_ang(i));    
  sigma2_Ibranch_I(i) = (Ibranch_mag(i)^2)*((cos(Ibranch_ang(i)))^2)*(sigma_thetaIbranch(i)^2) +...
                    (sigma_Ibranch(i)^2)*((sin(Ibranch_ang(i)))^2) + ...
                    Ibranch_mag(i)*sigma_Ibranch(i)*sigma_thetaIbranch(i)*sin(2*Ibranch_ang(i));     
end


%}

% Aggiunto il valore assoluto nell'ultimo termine come da email del
% 23/05/2020
sigma2_V_R = (V_mag.^2).*(sin(V_ang_rad).^2).*(sigma_thetaV.^2) + ...
             (sigma_V.^2).*(cos(V_ang_rad).^2) - ...
             V_mag.*sigma_V.*sigma_thetaV.*sin(2*V_ang_rad);
sigma2_V_I = (V_mag.^2).*(cos(V_ang_rad).^2).*(sigma_thetaV.^2) + ...
             (sigma_V.^2).*(sin(V_ang_rad).^2) + ...
             abs(V_mag.*sigma_V.*sigma_thetaV.*sin(2*V_ang_rad));

sigma2_Iinj_R = (Iinj_mag.^2).*(sin(Iinj_ang).^2).*(sigma_thetaIinj.^2) + ...
             (sigma_Iinj.^2).*(cos(Iinj_ang).^2) - ...
             Iinj_mag.*sigma_Iinj.*sigma_thetaIinj.*sin(2*Iinj_ang);
sigma2_Iinj_I = (Iinj_mag.^2).*(cos(Iinj_ang).^2).*(sigma_thetaIinj.^2) + ...
             (sigma_Iinj.^2).*(sin(Iinj_ang).^2) + ...
             abs(Iinj_mag.*sigma_Iinj.*sigma_thetaIinj.*sin(2*Iinj_ang));
         
sigma2_Ibranch_R = (Ibranch_mag.^2).*(sin(Ibranch_ang).^2).*(sigma_thetaIbranch.^2) + ...
             (sigma_Ibranch.^2).*(cos(Ibranch_ang).^2) - ...
             Ibranch_mag.*sigma_Ibranch.*sigma_thetaIbranch.*sin(2*Ibranch_ang);
sigma2_Ibranch_I = (Ibranch_mag.^2).*(cos(Ibranch_ang).^2).*(sigma_thetaIbranch.^2) + ...
             (sigma_Ibranch.^2).*(sin(Ibranch_ang).^2) + ...
             abs(Ibranch_mag.*sigma_Ibranch.*sigma_thetaIbranch.*sin(2*Ibranch_ang));
         
         
for i=1:nbranch
   if  sigma2_Ibranch_R(i)<=1e-10
       sigma2_Ibranch_R(i) =1e-10;       
   end
end
for i=1:nbranch
   if  sigma2_Ibranch_I(i)<=1e-10
       sigma2_Ibranch_I(i) =1e-10;       
   end
end

%{
    for i = 1:nbus
        if sigma2_V_R(i) > 10e-7
            sigma2_V_R(i) = sigma2_V_R(i);
        else
            sigma2_V_R(i) = 10e-7;
        end
        
        if sigma2_V_I(i) > 10e-7
            sigma2_V_I(i) = sigma2_V_I(i);
        else
            sigma2_V_I(i) = 10e-7;
        end
        
        if sigma2_I_R(i) > 10e-7
            sigma2_I_R(i) = sigma2_I_R(i);
        else
            sigma2_I_R(i) = 10e-7;
        end
        
        if sigma2_I_I(i) > 10e-7
            sigma2_I_I(i) = sigma2_I_I(i);
        else
            sigma2_I_I(i) = 10e-7;
        end
    end 
    %}
sigma2_Iinj_R = sigma2_Iinj_R*0 + 1e-7;
sigma2_Iinj_I = sigma2_Iinj_I*0 + 1e-7;
sigma2_Ibranch_R = sigma2_Ibranch_R*0 + 1e-7;
sigma2_Ibranch_I = sigma2_Ibranch_I*0 + 1e-7;
                 
% Deviazioni standard in coordinate rettangolari
std_V_R = sqrt(sigma2_V_R);
std_V_I = sqrt(sigma2_V_I);
std_Iinj_R = sqrt(sigma2_Iinj_R);
std_Iinj_I = sqrt(sigma2_Iinj_I);
std_Ibranch_R = sqrt(sigma2_Ibranch_R);
std_Ibranch_I = sqrt(sigma2_Ibranch_I);

% Matrice R    
R_V_R = sigma2_V_R;
R_V_I = sigma2_V_I;
R_Iinj_R = sigma2_Iinj_R;
R_Iinj_I = sigma2_Iinj_I;
R_Ibranch_R = sigma2_Ibranch_R;
R_Ibranch_I = sigma2_Ibranch_I;


% Matrice W
W_V_R = zeros(length(sigma2_V_R), 1);
W_V_I = zeros(length(sigma2_V_I), 1);
W_Iinj_R = zeros(length(sigma2_Iinj_R), 1);
W_Iinj_I = zeros(length(sigma2_Iinj_I), 1);
W_Ibranch_R = zeros(length(sigma2_Ibranch_R), 1);
W_Ibranch_I = zeros(length(sigma2_Ibranch_I), 1);
for kk = 1:length(sigma2_V_R)
    W_V_R(kk) = 1/sigma2_V_R(kk);
    W_V_I(kk) = 1/sigma2_V_I(kk);
    W_Iinj_R(kk) = 1/sigma2_Iinj_R(kk);
    W_Iinj_I(kk) = 1/sigma2_Iinj_I(kk);   
end
for kk = 1:length(sigma2_Ibranch_R)
    W_Ibranch_R(kk) = 1/sigma2_Ibranch_R(kk);
    W_Ibranch_I(kk) = 1/sigma2_Ibranch_I(kk);
end
    
% Varianza per misure virtuali Zero Injection
if ZIcount > 1
    min_sigma2_IR = min(sigma2_Iinj_R(setdiff(1:end,ZIbus)));
    R_ZI_IR = min_sigma2_IR*1e-2 * ones(numel(ZIbus),1);
    min_sigma2_II = min(sigma2_Iinj_I(setdiff(1:end,ZIbus)));
    R_ZI_II = min_sigma2_II*1e-2 * ones(numel(ZIbus),1);

    for i = 1:length(R_ZI_IR)
        W_ZI_IR(i,1) = 1/R_ZI_IR(i);
        W_ZI_II(i,1) = 1/R_ZI_II(i);
    end
else
    W_ZI_IR = [];
    W_ZI_II = [];
    R_ZI_IR = [];
    R_ZI_II = [];
end







end