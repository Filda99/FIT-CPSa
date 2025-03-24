clc; 
clear;

%{
Vysvetlivky:
Úhel náběhu (α): Úhel mezi podélnou osou letadla a směrem proudění vzduchu. Kladný úhel náběhu znamená, že nos letadla směřuje nahoru.
Úhel bočního skluzu (β): Úhel mezi podélnou osou letadla a směrem proudění vzduchu v horizontální rovině. Kladný úhel bočního skluzu znamená, že letadlo směřuje doprava.
Úhel klonění (φ): Úhel otočení letadla kolem podélné osy. Kladný úhel klonění znamená, že pravé křídlo směřuje dolů.
Úhel stoupání (θ): Úhel otočení letadla kolem příčné osy. Kladný úhel stoupání znamená, že nos letadla směřuje nahoru.
Úhel stáčení (ψ): Úhel otočení letadla kolem svislé osy. Kladný úhel stáčení znamená, že nos letadla směřuje doprava.

Rychlost klopení (p): Rychlost otáčení letadla kolem podélné osy.
Rychlost stoupání (q): Rychlost otáčení letadla kolem příčné osy.
Rychlost stáčení (r): Rychlost otáčení letadla kolem svislé osy.
Boční síla (Y): Síla působící na letadlo v bočním směru.
Klopivý moment (L): Moment, který otáčí letadlo kolem podélné osy.
Stáčivý moment (N): Moment, který otáčí letadlo kolem svislé osy.

Pr:
C_lbeta: Derivace klopivého momentu vzhledem k úhlu bočního skluzu. Vyjadřuje, jak moc se změní klopivý moment letadla, pokud se změní úhel bočního skluzu.
C_nbeta: Derivace stáčivého momentu vzhledem k úhlu bočního skluzu. Vyjadřuje, jak moc se změní stáčivý moment letadla, pokud se změní úhel bočního skluzu.
Vysvetleni:
    Derivace protože vyjadřuje míru změny jedné veličiny v závislosti na změně
    jiné veličiny. Říká nám, jak moc se změní klopivý moment letadla, pokud se změní úhel bočního skluzu o malou hodnotu.


%}
% =========================================================================
% Tabulka 4: Referenční geometrie a data o hmotnosti
S = 510.96;       % [m^2] Plocha křídla
c_bar = 8.32;     % [m] Střední aerodynamická tětiva
b = 59.74;        % [m] Rozpětí křídel
m = 288773.23;    % [kg] Hmotnost letadla
I_xx = 24675886.69; % [kg*m^2] Moment setrvačnosti kolem osy x
I_yy = 44877574.145; % [kg*m^2] Moment setrvačnosti kolem osy y
I_zz = 67384152.115; % [kg*m^2] Moment setrvačnosti kolem osy z
I_xz = 1315143.4115; % [kg*m^2] Produkt setrvačnosti
g = 9.81;           % Gravitacni konstanta

% =========================================================================
% Tabulka 5: Data letových podmínek a ustálené koeficienty stavu
h = 6096;         % [m] Nadmořská výška
Mach = 0.650;     % Machovo číslo
V_TAS_s = 205.13; % [m/s] Ustálená skutečná rychlost letu
VO_e = V_TAS_s;
q_bar = 13888;    % [N/m^2] Dynamický tlak
CG = 0.25;        % Poloha těžiště
alpha_s = 0.043633; % [rad] Úhel náběhu
C_L1 = 0.40;      % Koeficient vztlaku
C_D1 = 0.0250;    % Součinitel odporu
C_Tx1 = 0.0250;   % Tahový koeficient

% =========================================================================
% Tabulka 6: Podélné koeficienty a derivace stability
C_D0 = 0.0164;    % Součinitel odporu při nulovém vztlaku
C_Du = 0;         % Derivace odporu vzhledem k rychlosti
C_Dalpha = 0.20;  % Derivace odporu vzhledem k úhlu náběhu
C_Tx0 = -0.055;   % Tahový koeficient při nulové rychlosti
C_Txu = C_Tx1;      % TBD
C_L0 = 0.21;      % Součinitel vztlaku při nulovém úhlu náběhu
C_Lu = 0.13;      % Derivace vztlaku vzhledem k rychlosti
C_Lalpha = 4.4;   % Derivace vztlaku vzhledem k úhlu náběhu
C_Lalpha_dot = 7.0; % Derivace vztlaku vzhledem k časové změně úhlu náběhu
C_Lq = 6.6;       % Derivace vztlaku vzhledem k rychlosti stoupání
C_m0 = 0;         % Momentový koeficient při nulovém úhlu náběhu
C_m1 = C_m0;        % TBD
C_mu = 0.013;     % Derivace momentu vzhledem k rychlosti
C_malpha = -1.0;  % Derivace momentu vzhledem k úhlu náběhu
C_mdelta_alpha = -4.0; % Derivace momentu vzhledem k výchylce výškovky
C_mq = -20.5;      % Derivace momentu vzhledem k rychlosti stoupání
C_mTu = 0;          % TBD
C_mT1 = 0;          % TBD
C_mT_alpha = 0;          % TBD

% =========================================================================
% Tabulka 7: Derivace boční a směrové stability
C_lbeta = 0.0164; % Derivace klopivého momentu vzhledem k úhlu bočního skluzu
C_lp = 0;         % Derivace klopivého momentu vzhledem k rychlosti klopení
C_lr = 0.20;      % Derivace klopivého momentu vzhledem k rychlosti stáčení
C_ybeta = -0.055; % Derivace boční síly vzhledem k úhlu bočního skluzu
C_yp = 0.21;      % Derivace boční síly vzhledem k rychlosti klopení
C_yr = 0.13;      % Derivace boční síly vzhledem k rychlosti stáčení
C_nbeta = 4.4;    % Derivace stáčivého momentu vzhledem k úhlu bočního skluzu
C_nTbeta = 7.0;   % Derivace stáčivého momentu vzhledem k úhlu bočního skluzu vyvolaná tahem motoru.
C_np = 6.6;       % Derivace stáčivého momentu vzhledem k rychlosti klopení
C_nr = 0.0;       % Derivace stáčivého momentu vzhledem k rychlosti stáčení

% =========================================================================
% Tabulka 8: Derivace ovládání a momentu závěsu
C_ddelta_e = 0.0;    % Derivace odporu vzhledem k výchylce výškovky
C_ldelta_e = 0.32;   % Derivace vztlaku vzhledem k výchylce výškovky
C_mdelta_e= -1.30;   % Moment zavesu vyskovky v dusledku vychylky vyskovky
C_dih = 0;           % Indukovaný odpor výškovky
C_lih = -0.70;       % Indukovaný vztlak
C_mih = -2.7;        % Moment závěsu výškovky
C_ldelta_a = 0.013;  % Derivace klopivého momentu vzhledem k výchylce křidélek
C_ldelta_r = 0.008;  % Derivace klopivého momentu vzhledem k výchylce směrovky
C_ydelta_a = 0;      % Derivace boční síly vzhledem k výchylce křidélek
C_ydelta_r = 0.120;  % Derivace boční síly vzhledem k výchylce směrovky
C_ndelta_a = 0.0018; % Derivace stáčivého momentu vzhledem k výchylce křidélek
C_ndelta_r = -0.1;   % Derivace stáčivého momentu vzhledem k výchylce směrovky


% =========================================================================
% Tabulka 9: Podélné rozměrové derivace stability

% Derivace podélné síly vzhledem k rychlosti
X_u = (-q_bar * S * (C_Du + 2 * C_D1)) / (m * VO_e);       

% Derivace svislé síly vzhledem k rychlosti
Z_u = (-q_bar * S * (C_Lu + 2 * C_L1)) / (m * VO_e);        

% Derivace svislé síly vzhledem k výchylce výškovky
Z_delta_e = (-q_bar * S * C_ldelta_e) / m;                  

% Derivace podélné síly tahu vzhledem k rychlosti
X_Tu = (q_bar * S * (C_Txu + 2 * C_Tx1)) / (m * VO_e);      

% Derivace svislé síly vzhledem k úhlu náběhu
Z_alpha = (-q_bar * S * (C_Lalpha + C_D1)) / m;             

% Derivace momentu vzhledem k rychlosti
M_u = (q_bar * S * c_bar * (C_mu + 2 * C_m1)) / (I_yy * VO_e);          

% Derivace podélné síly vzhledem k úhlu náběhu
X_alpha = (-q_bar * S * (C_Dalpha - C_L1)) / m;             

% Derivace svislé síly vzhledem k časové změně úhlu náběhu
Z_alpha_dot = (-q_bar * S * c_bar * C_Lalpha_dot) / (2 * m * VO_e);     

% Derivace momentu tahu vzhledem k rychlosti
M_Tu = (q_bar * S * c_bar * (C_mTu + 2 * C_mT1)) / (I_yy * VO_e);       

% Derivace podélné síly vzhledem k výchylce výškovky
X_delta_e = (-q_bar * S * C_ddelta_e) / m;              

% Derivace svislé síly vzhledem k úhlové rychlosti příčného náklonu
Z_q = (-q_bar * S * c_bar * C_Lq) / (2 * m * VO_e);

% Derivace momentu vzhledem k úhlu náběhu
M_alpha = (q_bar * S * c_bar * C_malpha) / (I_yy);

% Derivace momentu vzhledem k časové změně úhlu náběhu
M_alpha_dot = (q_bar * S * c_bar^2 * C_Lalpha_dot) / (2 * I_yy * VO_e);

% Derivace momentu vzhledem k úhlové rychlosti příčného náklonu
M_q = (q_bar * S * c_bar^2 * C_mq) / (2 * I_yy * VO_e);

% Derivace momentu tahu vzhledem k úhlu náběhu
M_T_alpha = (q_bar * S * c_bar * C_mT_alpha) / I_yy;

% Derivace momentu vzhledem k výchylce výškovky
M_delta_e = (q_bar * S * c_bar * C_mdelta_e) / I_yy;


% =========================================================================
% Tabulka 10: Boční-směrové rozměrové derivace stability

% Derivace boční síly vzhledem k úhlu stranového skluzu (beta)
Y_beta = (q_bar * S * C_ybeta) / m;

% Derivace bočního momentu vzhledem k úhlu stranového skluzu (beta)
L_beta = (q_bar * S * b * C_lbeta) / I_xx;

% Derivace směrového momentu vzhledem k úhlu stranového skluzu (beta)
N_beta = (q_bar * S * b * C_nbeta) / I_zz;

% Derivace boční síly vzhledem k úhlové rychlosti klopení (p)
Y_p = (q_bar * S * b * C_yp) / (2 * m * VO_e);

% Derivace bočního momentu vzhledem k úhlové rychlosti klopení (p)
L_p = (q_bar * S * b^2 * C_lp) / (2 * I_xx * VO_e);

% Derivace bočního momentu tahu vzhledem k úhlu stranového skluzu (beta)
N_T_beta = (q_bar * S * b * C_nTbeta) / I_zz;

% Derivace boční síly vzhledem k úhlové rychlosti zatáčení (r)
Y_r = (q_bar * S * b * C_yr) / (2 * m * VO_e);

% Derivace bočního momentu vzhledem k úhlové rychlosti zatáčení (r)
L_r = (q_bar * S * b^2 * C_lr) / (2 * I_xx * VO_e);

% Derivace směrového momentu vzhledem k úhlové rychlosti klopení (p)
N_p = (q_bar * S * b^2 * C_np) / (2 * I_zz * VO_e);

% Derivace boční síly vzhledem k výchylce křidélek (delta_a)
Y_delta_a = (q_bar * S * C_ydelta_a) / m;

% Derivace bočního momentu vzhledem k výchylce křidélek (delta_a)
L_delta_a = (q_bar * S * b * C_ldelta_a) / I_xx;

% Derivace směrového momentu vzhledem k úhlové rychlosti zatáčení (r)
N_r = (q_bar * S * b^2 * C_nr) / (2 * I_zz * VO_e);

% Derivace boční síly vzhledem k výchylce směrovky (delta_r)
Y_delta_r = (q_bar * S * C_ydelta_r) / m;

% Derivace bočního momentu vzhledem k výchylce směrovky (delta_r)
L_delta_r = (q_bar * S * b * C_ldelta_r) / I_xx;

% Derivace směrového momentu vzhledem k výchylce křidélek (delta_a)
N_delta_a = (q_bar * S * b * C_ndelta_a) / I_zz;

% Derivace směrového momentu vzhledem k výchylce směrovky (delta_r)
N_delta_r = (q_bar * S * b * C_ndelta_r) / I_zz;

% Longitudinal equations in the matrixes

theta_s = 0; %TBD
phi_s   = 0;
psi_s   = 0;

P_s = 0.0;
Q_s = 0.0;
R_s = 0.0;

lon_A_pre = [ 
      (X_u + X_Tu),     X_alpha,                0,             -g*cos(theta_s);
      Z_u,            Z_alpha,                (Z_q + V_TAS_s), -g*sin(theta_s);
      (M_u + M_Tu),   (M_alpha + M_T_alpha),  M_q,             0;
      0,              0,                      1,               0 
];

lon_B_pre = [ 
    X_delta_e;
    Z_delta_e;
    M_delta_e;
    0 
];

lon_C = [ 
    1, 0,                   0, 0;
    0, (V_TAS_s - Z_alpha), 0, 0;
    0, -M_alpha_dot,        1, 0;
    0, 0,                   0, 1 
];

lon_A = lon_C' * lon_A_pre;
lon_B = lon_C' * lon_B_pre;


% Lateral directional state-space equations

lat_A_pre = [
    Y_beta,              Y_p, (Y_r - V_TAS_s), (g*cos(theta_s)), 0;
    L_beta,              L_p, L_r,             0,                0;
    (N_beta + N_T_beta), N_p, N_r,             0,                0;
    0,                   1,   0,               0,                0;
    0,                   0,   1,               0,                0
];

lat_B_pre = [
    Y_delta_a, Y_delta_r;
    L_delta_a, L_delta_r;
    N_delta_a, N_delta_r;
    0,         0;
    0,         0
];

lat_C = [
    V_TAS_s,  0,            0,            0, 0;
    0,        1,            (-I_xz/I_xx), 0, 0;
    0,        (-I_xz/I_zz), 1,            0, 0;
    0,        0,            0,            1, 0;
    0,        0,            0,            0, 1
];

lat_A = lat_C' * lat_A_pre;
lat_B = lat_C' * lat_B_pre;