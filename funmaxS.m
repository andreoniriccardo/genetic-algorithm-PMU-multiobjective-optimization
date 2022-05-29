function maxS = funmaxS(x)
global G Gbus_noZI Gbus_ZI g_line nnzGBbus_noZI nnzGBbus_ZI nnzgbline H_max_nominal R_max ZIbus
xGBbus_noZI_f = x(1:nnzGBbus_noZI);
xGBbus_ZI_f = x(nnzGBbus_noZI+1:nnzGBbus_noZI+nnzGBbus_ZI);
xgbline_f = x(nnzGBbus_noZI+nnzGBbus_ZI+1:nnzGBbus_noZI+nnzGBbus_ZI+nnzgbline);

sizex = nnzGBbus_noZI + nnzGBbus_ZI + nnzgbline;

M11_f = diag(ones(G,1));
M12_f = zeros(G,G);
M21_f = zeros(G,G);
M22_f = diag(ones(G,1));

[row_Gbus_noZI_f, col_Gbus_noZI_f] = find(Gbus_noZI);
MGbus_noZI_f = accumarray([row_Gbus_noZI_f(:),col_Gbus_noZI_f(:)],xGBbus_noZI_f(:));
MBbus_noZI_f = accumarray([row_Gbus_noZI_f(:),col_Gbus_noZI_f(:)],xGBbus_noZI_f(:));
MBbus_noZI_f(50,86)=1;
MBbus_noZI_f(49,87)=1;
MBbus_noZI_f(50,87)=1;

[row_Gbus_ZI_f, col_Gbus_ZI_f] = find(Gbus_ZI);
MGbus_ZI_f = accumarray([row_Gbus_ZI_f(:),col_Gbus_ZI_f(:)],xGBbus_ZI_f(:));
MBbus_ZI_f = accumarray([row_Gbus_ZI_f(:),col_Gbus_ZI_f(:)],xGBbus_ZI_f(:));
[abc, ncol_MGbus_ZI_f] = size(MGbus_ZI_f);
[abc, ncol_MBbus_ZI_f] = size(MBbus_ZI_f);
if ncol_MGbus_ZI_f < G
    diff = G-ncol_MGbus_ZI_f;
    MGbus_ZI_f(:,ncol_MGbus_ZI_f+1:G)=zeros(numel(ZIbus,diff));
end
if ncol_MBbus_ZI_f < G
    diff = G-ncol_MBbus_ZI_f;
    MBbus_ZI_f(:,ncol_MBbus_ZI_f+1:G)=zeros(numel(ZIbus,diff));
end


[row_gline_f, col_gline_f] = find(g_line);
Mgline_f = accumarray([row_gline_f(:),col_gline_f(:)],xgbline_f(:));
Mbline_f = accumarray([row_gline_f(:),col_gline_f(:)],xgbline_f(:));
Mbline_f(86,86) = 1;
Mbline_f(86,87) = 1;


Tol_f = [M11_f              M12_f;
       M21_f              M22_f;
       MGbus_noZI_f      -MBbus_noZI_f;
       MBbus_noZI_f       MGbus_noZI_f;
       MGbus_ZI_f        -MBbus_ZI_f;
       MBbus_ZI_f         MGbus_ZI_f;
       Mgline_f          -Mbline_f
       Mbline_f           Mgline_f];


H_tol = H_max_nominal.*Tol_f;
R_tilde = (1/0.0033^2)*R_max;
S = inv(H_tol'/R_tilde*H_tol);

logical = x<=0.9 | x>=1.1;
penalty = 1e10*ones(sizex,1);
%maxS = -max(max(S)) + penalty'*logical;
maxS = -max(max(S))
end