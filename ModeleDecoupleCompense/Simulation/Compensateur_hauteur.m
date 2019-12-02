%% Compensateur hauteur de la plaque
close all

% Specification
wg_des  = 185;
PM_des  = 25;
erp     = -0.0004;

figure; margin(FTplaque(3,3));  
%Les specs ne sont pas atteintes, donc compensation

% Ajout d'un compensateur d'avance de phase pour réglé le régime permanent
FTBO = FTplaque(3,3);
[num den] = tfdata(FTBO, 'v');
K_des = (1/abs(polyval(num, wg_des*1i)/polyval(den, wg_des*1i)));

figure
[GM,PM, wp, wg] = margin(FTBO*K_des);
deltaPhase = PM_des - PM;
alpha = (1-sind(deltaPhase))/(1+sind(deltaPhase));
T = 1/(wg_des*sqrt(alpha));
Ka = K_des/sqrt(alpha);
z = -1/T;
p = -1/(alpha*T); 
numGa = [1 -z];
denGa = [1 -p];
Ga = tf(Ka*numGa, denGa);

FTBF = FTBO * Ga;

%Cascader le compensateur AvPh avec FTplaque en z
FTplaque_z_comp = series(Ga_z, FTplaque(3,3));
[num_z_comp, den_z_comp] = tfdata(FTplaque_z_comp, 'v');

%Vérification des performances
disp('Performances de la FTplaque hauteur z compensé avec AvPh')
stepinfo(feedback(FTplaque_z_comp,1))
errRPech_z_comp = den_z_comp(end)/num_z_comp(end);

%Choisir entre le cas 1 ou le cas 2
%RePh (cas 2 oû l'erreur ech = -0.0004)
K_pos_z = 1/errRPech_z_comp;
K_pos_des_z = 1/errRPech_des_cas1_z;
K_des_z = K_pos_des_z/K_pos_z;
beta_z = abs(polyval(K_des_z*num_z_comp, wg_des_z*1i)/polyval(den_z_comp, wg_des_z*1i));
T_RePh_z = 10/wg_des_z;
z_RePh_z = -1/T_RePh_z;
p_RePh_z = -1/(beta_z*T_RePh_z);
Kr_z = K_des_z /beta_z;
numRePh_z = [T_RePh_z 1];
denRePh_z = [beta_z*T_RePh_z 1];
G_RePh_z = tf(Kr_z*beta_z*numRePh_z, denRePh_z);

% %PI (cas 2 oû l'erreur ech = 0)
% z_PI_z = -wg_des_z/10;
% numPI_z = [1 -z_PI_z];
% denPI_z = [1 0];
% Kp_z = polyval(conv(denPI_z, den_z_comp), wg_des_z*1i)/polyval(conv(numPI_z, num_z_comp), wg_des_z*1i);
% G_PI_z = tf(Kp2_z*numPI_z, denPI_z);

%Vérification des performances
%cas 1
FTplaque_z_comp21 = series(G_RePh_z, FTplaque_z_comp);
[num_z_comp21, den_z_comp21] = tfdata(FTplaque_z_comp21, 'v');
errRPech_z_comp21 = den_z_comp21(end)/num_z_comp21(end);

% %cas 2
% FTplaque_z_comp22 = series(G_PI_z, FTplaque_z_comp);
% [num_z_comp22, den_z_comp22] = tfdata(FTplaque_z_comp22, 'v');
% errRPech_z_comp22 = den_z_comp22(end)/num_z_comp22(end);