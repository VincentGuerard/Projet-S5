%% SPÉCIFICATIONS
%Inclinaison de la plaque (phi et theta)
Mp_des_angle = 0.05;
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
Ts_des_s = 3; %Entre 2 et 4
zeta_des_s = 0.9;   


%% CONCEPTION COMPENSATEURS
%*****Inclinaision de la plaque (phi et theta)*****%
phi_des_angle = atan(-pi/log(Mp_des_angle)); %rad
zeta_des_angle = cos(phi_des_angle);
wn1_des_angle = pi/(Tp_des_angle*sqrt(1-zeta_des_angle^2));
wn2_des_angle = (1+1.1*zeta_des_angle+1.4*zeta_des_angle^2)/Tr_des_angle;
wn3_des_angle = 4/(zeta_des_angle*Ts_des_angle);
wn_des_angle = max([wn1_des_angle wn2_des_angle wn3_des_angle]);
pole_des_angle = [-zeta_des_angle*wn_des_angle+wn_des_angle*sqrt(1-zeta_des_angle^2)*1i;
                  -zeta_des_angle*wn_des_angle-wn_des_angle*sqrt(1-zeta_des_angle^2)*1i];

figure; hold on; plot(pole_des_angle, 'p'); rlocus(FTplaque(1,1));  %Les pôles désirés sont placés à gauche                                                                   
title('Lieu des racines \phi / V_\phi et pôles désirés')            %des tracés dominant, donc un          
figure; hold on; plot(pole_des_angle, 'p'); rlocus(FTplaque(1,1));  %compensateur est nécessaire
title('Lieu des racines \theta / V_\theta et pôles désirés')

Ga_angle_phi = avPhase(pole_des_angle, FTplaque(1,1), 100*pi/180, 0);
Ga_angle_theta = avPhase(pole_des_angle, FTplaque(2,2), 100*pi/180, 0); 
FTplaque_phi_comp = series(FTplaque(1,1), Ga_angle_phi);
FTplaque_theta_comp = series(FTplaque(2,2), Ga_angle_theta);

%Affichage de l'atteinte des pôles désirés
figure; hold on; rlocus(FTplaque_phi_comp); plot(pole_des_angle, 'p'); hold off
title('Lieu des racines de \phi / V_\phi avec AvPh')
figure; hold on; rlocus(FTplaque_theta_comp); plot(pole_des_angle, 'p'); hold off
title('Lieu des racines de \theta / V_\theta avec AvPh')

%Vérification des performances

%*****Hauteur de la plaque*****%
figure; margin(FTplaque(3,3));  %Les specs ne sont pas atteintes, donc compensation
K_des_z = abs(polyval(denZ_VZ, wg_des_z*1i)/polyval(numZ_VZ, wg_des_z*1i));
figure; margin(FTplaque(3,3)*K_des_z);
[GM,PM, wp, wg] = margin(FTplaque(3,3)*K_des_z);


