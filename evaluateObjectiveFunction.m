function [f_eval, err] = pmu_evalZI_correnteLinee(chrom)

    global Incidence G E x_nom n_meas Gbus Bbus gen max_f2
    global R_V_R R_V_I R_Iinj_R R_Iinj_I W_V_R W_V_I W_Iinj_R W_Iinj_I V_R V_I Iinj_R Iinj_I
    global std_V_R std_V_I std_Iinj_R std_Iinj_I
    global ZIbus Gbus_ZI Bbus_ZI W_ZI_IR W_ZI_II R_ZI_IR R_ZI_II omega nbranch Branch_connectivity
    global C_line g_line b_line Ibranch_R Ibranch_I R_Ibranch_R R_Ibranch_I
    global W_Ibranch_R W_Ibranch_I std_Ibranch_R std_Ibranch_I cost n_sim_sens delta_H
    global nnzGBbus_noZI
    global M11_f3 M12_f3 M21_f3 M22_f3 Gbus_NoZI_f3 Bbus_NoZI_f3 Gbus_ZI_f3 Bbus_ZI_f3

    noZIbus = setdiff(1:G,ZIbus);
    % Number of bus voltage measurements
    tot_meas = nnz(chrom);
    % Number of injected current measurements
    tot_meas_noZI = nnz(chrom(setdiff(1:end,ZIbus)));
    % Number of line current measurements
    measured_branches = abs(C_line)*chrom';
    tot_meas_Ibranch = nnz(measured_branches);
    measured_branches_logical = 1*(measured_branches ~= 0);

    f1 = dot(cost, chrom');
    O = Incidence*chrom';

    % Initialization of measurement matrix
    V_meas_R = zeros(tot_meas,n_meas);
    V_meas_I = zeros(tot_meas,n_meas);
    Iinj_meas_R = zeros(tot_meas_noZI,n_meas);
    Iinj_meas_I = zeros(tot_meas_noZI,n_meas);
    Ibranch_meas_R = zeros(tot_meas_Ibranch,n_meas);
    Ibranch_meas_I = zeros(tot_meas_Ibranch,n_meas);

    % Initialization of matrix H
    M11 = zeros(tot_meas,G);
    M12 = zeros(tot_meas,G);
    M21 = zeros(tot_meas,G);
    M22 = zeros(tot_meas,G);
    G_new = zeros(tot_meas_noZI,G);
    B_new = zeros(tot_meas_noZI,G);
    g_new = zeros(tot_meas_Ibranch,G);
    b_new = zeros(tot_meas_Ibranch,G);

    % Initialization of covariance matrix
    W_V_R_new = zeros(tot_meas,1);
    W_V_I_new = zeros(tot_meas,1);
    W_Iinj_R_new = zeros(tot_meas_noZI,1);
    W_Iinj_I_new = zeros(tot_meas_noZI,1);
    W_Ibranch_R_new = zeros(tot_meas_Ibranch,1);
    W_Ibranch_I_new = zeros(tot_meas_Ibranch,1);
    R_V_R_new = zeros(tot_meas,1);
    R_V_I_new = zeros(tot_meas,1);
    R_Iinj_R_new = zeros(tot_meas_noZI,1);
    R_Iinj_I_new = zeros(tot_meas_noZI,1);
    R_Ibranch_R_new = zeros(tot_meas_Ibranch,1);
    R_Ibranch_I_new = zeros(tot_meas_Ibranch,1);

    for meas = 1:n_meas
        count = 1;
        count_noZI = 1;
        countIbranch = 1;

        % Simulation of measurements
        for i = 1:G
            if ~ismember(i, ZIbus)
                for k = 1:chrom(i)
                    % 0.0001 to avoid singularity
                    V_meas_R(count,meas) = normrnd(V_R(i), abs(std_V_R(i)));
                    V_meas_I(count,meas) = normrnd(V_I(i), abs( ...
                                                           std_V_I(i)));
                    Iinj_meas_R(count_noZI,meas) = normrnd(Iinj_R(i), abs( ...
                                                           std_Iinj_R(i)));
                    Iinj_meas_I(count_noZI,meas) = normrnd(Iinj_I(i), abs( + ...
                                                           std_Iinj_I(i)));

                    W_V_R_new(count,1) = W_V_R(i);
                    W_V_I_new(count,1) = W_V_I(i);
                    W_Iinj_R_new(count_noZI,1) = W_Iinj_R(i);
                    W_Iinj_I_new(count_noZI,1) = W_Iinj_I(i);
                    R_V_R_new(count,1) = R_V_R(i);
                    R_V_I_new(count,1) = R_V_I(i);
                    R_Iinj_R_new(count_noZI,1) = R_Iinj_R(i);
                    R_Iinj_I_new(count_noZI,1) = R_Iinj_I(i);

                    M11(count,i) = 1;
                    M22(count,i) = 1;
                    G_new(count_noZI,:) = Gbus(i,:);
                    B_new(count_noZI,:) = Bbus(i,:);

                    count = count + 1;
                    count_noZI = count_noZI + 1;
                end

            else
                for k = 1:chrom(i)
                    V_meas_R(count,meas) = normrnd(V_R(i), abs(0.000000001 + ...
                                                           std_V_R(i)));
                    V_meas_I(count,meas) = normrnd(V_I(i), abs(0.000000001 + ...
                                                           std_V_I(i)));

                    W_V_R_new(count,1) = W_V_R(i);
                    W_V_I_new(count,1) = W_V_I(i);
                    R_V_R_new(count,1) = R_V_R(i);
                    R_V_I_new(count,1) = R_V_I(i);

                    M11(count,i) = 1;
                    M22(count,i) = 1;

                    count = count + 1;
                end
            end
        end

        for i = 1:nbranch
            for k = 1:measured_branches_logical(i)
                Ibranch_meas_R(countIbranch,meas) = normrnd(Ibranch_R(i), abs(0.000000001 + ...
                                                           std_Ibranch_R(i)));
                Ibranch_meas_I(countIbranch,meas) = normrnd(Ibranch_I(i), abs(0.000000001 + ...
                                                           std_Ibranch_I(i)));

                W_Ibranch_R_new(countIbranch,1) = W_Ibranch_R(i);
                W_Ibranch_I_new(countIbranch,1) = W_Ibranch_I(i);
                R_Ibranch_R_new(countIbranch,1) = R_Ibranch_R(i);
                R_Ibranch_I_new(countIbranch,1) = R_Ibranch_I(i);

                g_new(countIbranch,:) = g_line(i,:);
                b_new(countIbranch,:) = b_line(i,:);

                countIbranch = countIbranch + 1;
            end
        end

    end

    % Add the ZI bus equations
    if numel(ZIbus) > 0
        I_ZI_R = zeros(numel(ZIbus),n_meas);
        I_ZI_I = zeros(numel(ZIbus),n_meas);
        Z = [V_meas_R; V_meas_I; Iinj_meas_R; Iinj_meas_I; I_ZI_R; I_ZI_I; Ibranch_meas_R; Ibranch_meas_I];

        H = [M11 M12;
             M21 M22;
             G_new    -B_new;
             B_new     G_new;
             Gbus_ZI  -Bbus_ZI;
             Bbus_ZI   Gbus_ZI;
             g_new    -b_new;
             b_new     g_new];

        W = diag([W_V_R_new; W_V_I_new; W_Iinj_R_new; W_Iinj_I_new; W_ZI_IR; W_ZI_II; W_Ibranch_R_new; W_Ibranch_I_new]);
    else
        Z = [V_meas_R; V_meas_I; Iinj_meas_R; Iinj_meas_I; Ibranch_meas_R; Ibranch_meas_I];
        H = [M11 M12;
             M21 M22;
             G_new    -B_new;
             B_new     G_new;
             g_new    -b_new;
             b_new     g_new];
        W = diag([W_V_R_new; W_V_I_new; W_Iinj_R_new; W_Iinj_I_new; W_Ibranch_R_new; W_Ibranch_I_new]);
    end

    Temp = H'*W*H; % Temporary matrix
    % Inverse of H'*W*H
    invTemp = inv(Temp);
    if sum(sum(isinf(invTemp)))==0 && sum(sum(isnan(invTemp))) == 0
    x_est = zeros(2*G, n_meas);
    est_error = zeros(2*G, n_meas);
    for ii = 1:n_meas
        x_est(:,ii) = (invTemp)*H'*(W)*Z(:,ii);
        est_error(:,ii) = x_nom - x_est(:,ii);
    end

    % Complex matrix of measurement errors
    est_error_new = zeros(G, n_meas);
    for i=1:G
        est_error_new(i,:) = est_error(i,:) + 1i*est_error(i+G,:);
    end
    % Total covariance matrix
    C = (1/n_meas)*(est_error_new*est_error_new');

   .% Measurement uncertainty
    f2 = max(sqrt(abs(eig(C))));

    %% f_3: Sensitivity
    size_noZI = G - numel(ZIbus);

    VR_index = find(chrom);
    VR_index_trasl = VR_index;
    VI_index = find(chrom);
    VI_index_trasl = G+VI_index;

    IinjR_index = setdiff(VR_index,ZIbus);
    IinjR_index_trasl = (2*G) + IinjR_index;
    IinjI_index = setdiff(VI_index,ZIbus);
    IinjI_index_trasl = (2*G) + size_noZI + IinjI_index;

    ZIR_index = 1:numel(ZIbus);
    ZIR_index_trasl = (2*G) + 2*size_noZI + ZIR_index;
    ZII_index = 1:numel(ZIbus);
    ZII_index_trasl = ((2*G) + 2*size_noZI + numel(ZIbus)) + ZII_index;

    IbranchR_index = find(measured_branches_logical)';
    IbranchR_index_trasl = ((2*G) + 2*size_noZI + 2*numel(ZIbus)) + IbranchR_index;
    IbranchI_index = find(measured_branches_logical)';
    IbranchI_index_trasl = ((2*G) + 2*size_noZI + 2*numel(ZIbus) + nbranch) + IbranchI_index;

    index_union = [VR_index_trasl VI_index_trasl IinjR_index_trasl IinjI_index_trasl ZIR_index_trasl ZII_index_trasl IbranchR_index_trasl IbranchI_index_trasl];

    % Covariance matrix
    R = diag([R_V_R_new; R_V_I_new; R_Iinj_R_new; R_Iinj_I_new; R_ZI_IR; R_ZI_II; R_Ibranch_R_new; R_Ibranch_I_new]);
    R_tilde = (1/0.0033^2)*R;

    H_f3 = delta_H(index_union,:);

    S_f3 = inv(H_f3'/R_tilde*H_f3);

    f3 = max(max(S_f3));

    %% Evaluates functions:
    % f1 Number of PMUs
    % f2 Measurement uncertainty
    % f3 Sensitivity
    f_eval = [f1 f2 f3];
    %% Observability error

    if f2 <= max_f2(gen)
        u = O>=1;

        v = zeros(G,1);
        for j = 1:numel(ZIbus)
            if Incidence(ZIbus(j),:)*u >= (omega(j)-1)
                v(ZIbus(j)) = 1;
            end
        end

        c1 = ones(G,1) - Incidence*u - Incidence*v;
        err1 = (c1>0).*c1;

        % PMU contingencies constraints
        c2 = zeros(G*nbranch,1);

        for ind_i = 1:nbranch
            bus1 = Branch_connectivity(ind_i,1);
            bus2 = Branch_connectivity(ind_i,2);
            Al = Incidence;
            Al(bus1,bus2) = 0;
            Al(bus2,bus1) = 0;
            Ol = Al*chrom';

            ul = Ol>=1;
            omegal = zeros(numel(ZIbus),1);
            for ind_j1 = 1:numel(ZIbus)
               omegal(ind_j1) = nnz(Al(ZIbus(ind_j1),:));
            end

            vl = zeros(G,1);
            for ind_j2 = 1:numel(ZIbus)
                if Al(ZIbus(ind_j2),:)*ul >= (omegal(ind_j2)-1)
                    vl(ZIbus(ind_j2)) = 1;
                end
            end
            % Line contingencies constraints evaluation
            c2(1+G*(ind_i-1): ind_i*G,1) = ones(G,1) - Al*ul - Incidence*vl;
        end
        err2 = (c2>0).*c2;

        % PMU contingencies constraints

        for i = 1:G

            Ap = Incidence;
            Ap(:,i) = zeros(G,1);

            Op = Ap*chrom';

            up = Op>=1;

            vp = zeros(G,1);

            for j = 1:numel(ZIbus)
                if Incidence(ZIbus(j),:)*up >= (omega(j)-1)
                    vp(ZIbus(j)) = 1;
                end
            end
            % PMU contingencies constraints evaluation
            c3(1+G*(i-1): i*G , 1) = ones(G,1) - Ap*up - Incidence*vp;
        end
        err3 = (c3>0).*c3;

        err = [err1', err2', err3'];
    else
        err = ones(1,E);
    end
    else
        f1 = Inf;
        f2 = Inf;
        f3 = Inf;
        f_eval = [f1 f2 f3];
        err = ones(1,E);
    end

end
