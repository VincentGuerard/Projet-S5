%% Compensateur hauteur Z
close all

% Specifications
PM_des_z            = 25;        %degré
wg_des_z            = 185;       %rad/sec
errRPech_des_cas1_z = -0.0004;   %Pour une consigne 0.01
errRPech_des_cas2_z = 0;
% Vérification des critères de performance
figure
margin(FTplaque(3,3));
title('Vérification des critères de performance FT original')
grid on

% Ajout d'un avance de phase Bode
K_des_z = (1./abs(polyval(numZ_VZ, wg_des_z*j)/polyval(denZ_VZ, wg_des_z*j)));
% figure
% margin(FTplaque(3,3)*K_des_z);
[GM,PM, wp, wg] = margin(FTplaque(3,3)*K_des_z);
deltaPhase_z = PM_des_z - PM;
alpha_z = (1-sind(deltaPhase_z))/(1+sind(deltaPhase_z));
T_z = 1/(wg_des_z*sqrt(alpha_z));
Ka_z = K_des_z/sqrt(alpha_z);
z_z = -1/T_z;
p_z = -1/(alpha_z*T_z); 
numGa_z = [1 -z_z];
denGa_z = [1 -p_z];
Ga_z = tf(Ka_z*numGa_z, denGa_z);

%Cascader le compensateur AvPh avec FTplaque en z
FTplaque_z_comp = series(Ga_z, FTplaque(3,3))
[num_z_comp, den_z_comp] = tfdata(FTplaque_z_comp, 'v');

% Ajout d'un 2e avance de phase Bode
K_des_z2 = (1./abs(polyval(num_z_comp, wg_des_z*j)/polyval(den_z_comp, wg_des_z*j)));
[GM2,PM2, wp2, wg2] = margin(FTplaque_z_comp*K_des_z2);
deltaPhase_z2 = PM_des_z - PM2 + 6; % Ajout de 6
alpha_z2 = (1-sind(deltaPhase_z2))/(1+sind(deltaPhase_z2));
T_z2 = 1/(wg_des_z*sqrt(alpha_z2));
Ka_z2 = K_des_z2/sqrt(alpha_z2);
z_z2 = -1/T_z2;
p_z2 = -1/(alpha_z2*T_z2); 
numGa_z2 = [1 -z_z2];
denGa_z2 = [1 -p_z2];
Ga_z2 = tf(Ka_z2*numGa_z2, denGa_z2);

%Cascader le 2e compensateur AvPh avec FTplaque en z
FTplaque_z_comp = series(Ga_z2, FTplaque_z_comp);
[num_z_comp, den_z_comp] = tfdata(FTplaque_z_comp, 'v');

%Vérification des performances
disp('Performances de la FTplaque hauteur z compensé avec AvPh')
stepinfo(feedback(FTplaque_z_comp,1))
errRPech_z_comp = den_z_comp(end)/num_z_comp(end);
[GM,PM, wp, wg] = margin(FTplaque_z_comp);
figure(11)
margin(FTplaque_z_comp)

%Choisir entre le cas 1 ou le cas 2
%RePh (cas 1 oû l'erreur ech = -0.0004)
K_pos_z = 1/errRPech_z_comp;
K_pos_des_z = 1/errRPech_des_cas1_z;
K_des_z = K_pos_des_z/K_pos_z;
beta_z = abs(polyval(K_des_z*num_z_comp, wg_des_z*1i)/polyval(den_z_comp, wg_des_z*1i));
T_RePh_z = 11/wg_des_z;
z_RePh_z = -1/T_RePh_z;
p_RePh_z = -1/(beta_z*T_RePh_z);
Kr_z = K_des_z /beta_z;
numRePh_z = [1 -z_RePh_z];
denRePh_z = [1 -p_RePh_z];
G_RePh_z = tf(Kr_z*numRePh_z, denRePh_z);
%PI (cas 2 oû l'erreur ech = 0)
z_PI_z = -wg_des_z/10;
numPI_z = [1 -z_PI_z];
denPI_z = [1 0];
Kp_z = abs(polyval(conv(denPI_z, den_z_comp), wg_des_z*1i)/polyval(conv(numPI_z, num_z_comp), wg_des_z*1i));
G_PI_z = tf(Kp_z*numPI_z, denPI_z);

%Vérification des performances
%cas 1
FTplaque_z_comp21 = series(G_RePh_z, FTplaque_z_comp);
[num_z_comp21, den_z_comp21] = tfdata(FTplaque_z_comp21, 'v');
errRPech_z_comp21 = den_z_comp21(end)/num_z_comp21(end);
%cas 2
FTplaque_z_comp22 = series(G_PI_z, FTplaque_z_comp);
[num_z_comp22, den_z_comp22] = tfdata(FTplaque_z_comp22, 'v');
errRPech_z_comp22 = den_z_comp22(end)/num_z_comp22(end);
figure(12)
margin(FTplaque_z_comp22)

FT_comp_z = Ga_z*Ga_z2*G_PI_z;
[numCompz, denCompz] = tfdata(FT_comp_z, 'v');
