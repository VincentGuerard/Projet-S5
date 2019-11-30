%% SPÉCIFICATIONS
%Inclinaison de la plaque (phi et theta)
Mp_des_angle = 0.0005;
Ts_des_angle = 0.03;
Tp_des_angle = 0.025;
Tr_des_angle = 0.02;
errRPech_des_angle = 0;

%Hauteur de la plaque (z)
PM_des_z = 25;    %degré
wg_des_z = 185;   %rad/sec
errRPech_des_cas1_z = -0.0004;    %Pour une consigne 0.01
errRPech_des_cas2_z = 0;

%Position de la sphère (Xs et Ys)
Ts_des_s    = 3; %Entre 2 et 4
zeta_des_s  = 0.9;   

%% CONCEPTION COMPENSATEURS
close all
clc
%**************************************************%
%********* Inclinaision de la plaque **************%
%**************************************************%
%Avance de phase
phi_des_angle   = atan(-pi/log(Mp_des_angle)); %rad
zeta_des_angle  = cos(phi_des_angle);
wn1_des_angle   = pi/(Tp_des_angle*sqrt(1-zeta_des_angle^2));
wn2_des_angle   = (1+1.1*zeta_des_angle+1.4*zeta_des_angle^2)/Tr_des_angle;
wn3_des_angle   = 4/(zeta_des_angle*Ts_des_angle);
wn_des_angle    = max([wn1_des_angle wn2_des_angle wn3_des_angle]);
pole_des_angle  = [-zeta_des_angle*wn_des_angle+wn_des_angle*sqrt(1-zeta_des_angle^2)*1i;
                   -zeta_des_angle*wn_des_angle-wn_des_angle*sqrt(1-zeta_des_angle^2)*1i];

figure; hold on; plot(pole_des_angle, 'p'); rlocus(FTplaque(1,1));  %Les pôles désirés sont placés à gauche                                                                   
title('Lieu des racines \phi / V_\phi et pôles désirés')            %des tracés dominant, donc un compensateur est necessaire

Ga_angle                         = avPhase(pole_des_angle, FTplaque(1,1), (15*pi/180), 2);
FTplaque_angles_comp             = series(FTplaque(1,1),Ga_angle);
[num_angle_comp, den_angle_comp] = tfdata(FTplaque_angles_comp, 'v');

Ga_angles = avPhase(pole_des_angle, FTplaque(1,1), (15/360)*2*pi, 2);
FTplaque_angles_comp = FTplaque(1,1) * 3 * Ga_angles;
info = stepinfo(feedback(FTplaque_angles_comp,1));

%Affichage de l'atteinte des pôles désirés
figure;
rlocus(FTplaque_angles_comp); hold on;
plot(pole_des_angle, 'p');
p = rlocus(FTplaque_angles_comp, 1);
plot(real(p), imag(p), 's'); hold off
xlim([-600 200]); ylim([-600 600])
title('Lieu des racines de \phi / V_\phi avec AvPh')

%Vérification des performances
disp('Performances FTplaque des angles compensée AvPh:')
stepinfo(feedback(FTplaque_angles_comp, 1))
t = 0:0.00001:0.25;
u = ones(length(t),1);
figure
lsim(feedback(FTplaque_angles_comp,1),u,t)

%PI
z_PI_angle = real(pole_des_angle(1))/5;
num_PI_angle = [1 -z_PI_angle];
den_PI_angle = [1 0];
Kp_angle = abs(polyval(conv(den_PI_angle, den_angle_comp), pole_des_angle(1))/polyval(conv(num_PI_angle, num_angle_comp), pole_des_angle(1)));
G_PI_angle = tf(1.5*Kp_angle*num_PI_angle, den_PI_angle);

%Cascader le PI à la FT compensée avec AvPh
FTplaque_angle_comp2 = series(FTplaque_angles_comp, G_PI_angle);
[num_angle_comp2, den_angle_comp2] = tfdata(FTplaque_angle_comp2, 'v');

% Affichage de l'atteinte des poles desirer
figure;
rlocus(FTplaque_angle_comp2); hold on;
plot(pole_des_angle, 'p');
p = rlocus(FTplaque_angle_comp2, 1);
plot(real(p), imag(p), 's'); hold off
xlim([-600 200]); ylim([-600 600])
title('Lieu des racines de \phi / V_\phi avec AvPh+PI')

%Vérification des performances
errRPech_angle = den_angle_comp2(end)/num_angle_comp2(end-1)
disp('Performances FTplaque des angles compensée AvPh et PI:')
stepinfo(feedback(FTplaque_angle_comp2, 1))

figure
lsim(feedback(FTplaque_angle_comp2,1),u,t);

%**************************************************%
%************ Hauteur de la plaque (z) ************%
%**************************************************%
% figure; margin(FTplaque(3,3));  %Les specs ne sont pas atteintes, donc compensation
% K_des_z = (1/abs(polyval(numZ_VZ, wg_des_z*1i)/polyval(denZ_VZ, wg_des_z*1i)));
% figure; margin(FTplaque(3,3)*K_des_z);
% [GM,PM, wp, wg] = margin(FTplaque(3,3)*K_des_z);
% deltaPhase_z = PM_des_z - PM;
% alpha_z = (1-sind(deltaPhase_z))/(1+sind(deltaPhase_z));
% T_z = 1/(wg_des_z*sqrt(alpha_z));
% Ka_z = K_des_z/sqrt(alpha_z);
% z_z = -1/T_z;
% p_z = -1/(alpha_z*T_z); 
% numGa_z = [1 -z_z];
% denGa_z = [1 -p_z];
% Ga_z = tf(Ka_z*numGa_z, denGa_z);
% 
% %Cascader le compensateur AvPh avec FTplaque en z
% FTplaque_z_comp = series(Ga_z, FTplaque(3,3));
% [num_z_comp, den_z_comp] = tfdata(FTplaque_z_comp, 'v');
% 
% %Vérification des performances
% disp('Performances de la FTplaque hauteur z compensé avec AvPh')
% stepinfo(feedback(FTplaque_z_comp,1))
% errRPech_z_comp = den_z_comp(end)/num_z_comp(end);
% 
% %Choisir entre le cas 1 ou le cas 2
% %RePh (cas 2 oû l'erreur ech = -0.0004)
% K_pos_z = 1/errRPech_z_comp;
% K_pos_des_z = 1/errRPech_des_cas1_z;
% K_des_z = K_pos_des_z/K_pos_z;
% beta_z = abs(polyval(K_des_z*num_z_comp, wg_des_z*1i)/polyval(den_z_comp, wg_des_z*1i));
% T_RePh_z = 10/wg_des_z;
% z_RePh_z = -1/T_RePh_z;
% p_RePh_z = -1/(beta_z*T_RePh_z);
% Kr_z = K_des_z /beta_z;
% numRePh_z = [T_RePh_z 1];
% denRePh_z = [beta_z*T_RePh_z 1];
% G_RePh_z = tf(Kr_z*beta_z*numRePh_z, denRePh_z);
% 
% % %PI (cas 2 oû l'erreur ech = 0)
% % z_PI_z = -wg_des_z/10;
% % numPI_z = [1 -z_PI_z];
% % denPI_z = [1 0];
% % Kp_z = polyval(conv(denPI_z, den_z_comp), wg_des_z*1i)/polyval(conv(numPI_z, num_z_comp), wg_des_z*1i);
% % G_PI_z = tf(Kp2_z*numPI_z, denPI_z);
% 
% %Vérification des performances
% %cas 1
% FTplaque_z_comp21 = series(G_RePh_z, FTplaque_z_comp);
% [num_z_comp21, den_z_comp21] = tfdata(FTplaque_z_comp21, 'v');
% errRPech_z_comp21 = den_z_comp21(end)/num_z_comp21(end);
% 
% % %cas 2
% % FTplaque_z_comp22 = series(G_PI_z, FTplaque_z_comp);
% % [num_z_comp22, den_z_comp22] = tfdata(FTplaque_z_comp22, 'v');
% % errRPech_z_comp22 = den_z_comp22(end)/num_z_comp22(end);
% 
% 
% 
% %**************************************************%
% %**************** Sphère (xs et ys) ***************%
% %**************************************************%