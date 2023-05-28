% Playin PMU data to Centralia Generator with subtransient model 
% May 2023

% *** initial v_mag v_angle p and q will be replaced***
bus = [ ...
% num volt               angle                p_gen                 q_gen                  p_load q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
  1   1.079696495811734  -25.797498136638520  -0.748043509171887    -0.172677797461532     0      0      0       0       1    100   -100  20      1.5   0.5;
  2   1.117753961796009  -19.428399188282760  0.750140458834272     0.263562087945019      0      0      0       0       2    100   -100  20      1.5   0.5];


line = [ ...
% bus bus r        x       y    tapratio tapphase tapmax tapmin tapsize
  1   2   0.00348  0.179742  0.0  1.0      0        0      0      0];%xformer

%Both generator parameters are from the example machine in chap 4 of
%Anderson and Fouad with a salient pole assumption.
%This model is a sub-transient model.
%From eqn (4.289), T"_qo=T'_qo if it is a 2-axis transient model.
mac_con = [ ...
% num bus base  NA   r_a    NA  x'_d    NA    Td     Tv       NA     NA    NA    NA    NA     H = 0  NA  NA bus
   1   1  100   0.0  0.0    0   1e-9    0     0.0009 0.0009   0      0     0.0   0.0   0.0    0      0   0   1   0      0;   % don't change (playin machine)
% num bus base  x_l  r_a    x_d  x'_d   x"_d  T'_do  T"_do    x_q    x'_q  x"_q  T'_qo T"_qo  H      d_0 d_1 bus S1     S12 
   2   2  100   0.12 0.0019 2.05 0.42   0.24  5.8    0.015    1.95   0.65  0.24  0.6   0.03   3.25   0   0   2   0.125  0.33]; % what about r and x comp (0,0.063)
  
%Exciter
%From p. 1137 of Kundur
% exc_con = [...
% %  type mach  T_R   K_A   T_A   T_B   T_C   Vmax  Vmin
%     0    2    0.02  200   0.02  0.02  0.32  10    -10]; %Fast static exciter, with TGR
% exc_con = [
% % Type Mac TR    KVP TA    KVI          KC   VRmax VRmin KE   TE   E1   SE1    E2   SE2   KVD   TVD   KD   Vimax  limflg
%   4    2   0.02  200 0.02  0.00000000   0.15 10    -10   1.0  1.5  6.5  0.3    9    3     60    0.02  0.45 999    0.0;];

%PSS
%Designed using the "large-inertia, infinite-bus" method in Roger's book and
%Kundur's book.
%type gen# K  Tw T1     T2     T3     T4     max  min
% pss_con = [ ...
%   1   2    2  5  0.4773 0.1675 0.4773 0.1675 0.05 -0.05];

% tg_con = [...
% %     mach wf 1/R     Tmax     Ts        Tc     T3     T4     T5
% %     1   2   1   20.0    1.00   0.40      45.0   5      -1     0.5];%hydro, twice as slow as steam
%    1   2   1   20.0    1.00   0.04      0.2    0.0    1.5    5.0];%steam, works good

% non-conforming load
% col 1       bus number
% col 2       fraction const active power load
% col 3       fraction const reactive power load
% col 4       fraction const active current load
% col 5       fraction const reactive current load
% load_con = [...
% %bus Pcont Qconst P_Icont Q_Iconst
%   8   1     1     0       0;
%   9   0     0     0       0];

%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  
%                  - 0 three phase
%                  - 1 line to ground
%                  - 2 line-to-line to ground
%                  - 3 line-to-line
%                  - 4 loss of line with no fault
%                  - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)

%% Auto replace this?
sw_con = [...
0         0    0    0    0    0    1/600;%sets intitial time step
28.8       1    2    0    0    6    1/600;
28.9       0    0    0    0    0    1/600;
29         0    0    0    0    0    1/600]; % end simulation

playin_con = {
    1 '/Users/ron/Documents/Ron/PhD Research/Code/pst_val/data/ex_playin_BPA_data/d_centralia_play_1223239to1225021.mat';
    0 ''};
