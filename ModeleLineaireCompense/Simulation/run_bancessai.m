close all
clear all
clc

% Position à l'équilibre de la sphère (pour tests statiques)
sig = 1.0;         % Présence (1) ou non (0) de la sphère
xSeq = 0.000;      % Position x de la sphère à l'équilibre en metres
ySeq = 0.000;      % Position y de la sphère à l'équilibre en metres

%Point d'opération choisi pour la plaque
Axeq = 0;               %en degres
Ayeq = 0;               %en degres
Pzeq = .015;            %en metres

%Exemple de trajectoire
t_des     = [0:1:8]'*5;

x_des     = [t_des, [0 0 0.5 1  0 -1 0 1 0]'*0.05];
y_des     = [t_des, [0 0 0 0 -1  0 1 0 0]'*0.05];
z_des     = [t_des, [1 1 1 1  1  1 1 1 1]'*.015];
tfin = 50;

%initialisation
bancEssaiConstantes

%bancessai_ini  %faites tous vos calculs de modele ici
run('../../Identification/Projet_s5.m');close all;
run('../Projet_S5_lin.m');

%Calcul des compensateurs
%iniCTL_ver4    %Calculez vos compensateurs ici
numCompAngles = 1.0e09*[0.000230814504577 0.026170635469275 0.897433010674818 9.699372248008062];
denCompAngles = 1.0e07*[0.000000100000000 0.000823402569667 1.694979479335234 0];
numCompz = 1.0e11*[0.000031609546854 0.004304100455063 0.142047621936623 1.354942422178642];
denCompz = 1.0e05*[0.000010000000000 0.017380342819025 5.055406372627656 0];

%simulation
open_system('DYNctl_ver4_etud_obfusc')
set_param('DYNctl_ver4_etud_obfusc','AlgebraicLoopSolver','TrustRegion')
sim('DYNctl_ver4_etud_obfusc')


%% affichage
% Comparaison de la position et de la vitesse de la bille
figure
subplot(211); 
hold on;
plot(tsim, ynonlineaire(:,7), 'b--', tsim, ynonlineaire(:,8), 'r--', 'LineWidth', 2); %Positions Non Linéaire
plot(tsim, ynonlineaire(:,9), 'g--', tsim, ynonlineaire(:,10), 'y--', 'LineWidth', 2); %Vitesses Non Linéaire
plot(tsim, yLineaire.Data(:,4), 'b', tsim, yLineaire.Data(:,5), 'r'); %Position Linéaire
plot(tsim, yLineaire.Data(:,6), 'g', tsim, yLineaire.Data(:,7), 'y'); %Vitesse Linéaire
legend('Px (NL)', 'Py (NL)', 'Vx (NL)', 'Vy (NL)', 'Px (L)', 'Py (L)', 'Vx (L)', 'Vy (L)');grid on; grid minor;
ylim([-0.06 0.06]);
title('Position et vitesse Non Linéaire (pointillé) vs Linéaire (ligne pleine)');

subplot(212); 
plot(tsim, yLineaire.Data(:,4) - ynonlineaire(:,7), 'b'); hold on; %Position x
plot(tsim, yLineaire.Data(:,5) - ynonlineaire(:,8), 'r'); %Position y
plot(tsim, yLineaire.Data(:,6) - ynonlineaire(:,9), 'g'); %Vitesse x
plot(tsim, yLineaire.Data(:,7) - ynonlineaire(:,10), 'y'); %Vitesse y
legend('\Delta Px', '\Delta Py', '\Delta Vx', '\Delta Vy');grid on; grid minor;
ylim([-0.02 0.02])
title('Différence entre Non Linéaire et Linéaire');
hold off;

% Comparaison des distances des capteurs
figure
subplot(211); 
hold on;title('Distance ZD ZE ZF Non Linéaire (Pointillé) vs Linéaire (Ligne pleine)');
plot(tsim, ynonlineaire(:,17), 'b--', 'LineWidth', 2);
plot(tsim, ynonlineaire(:,18), 'r--', 'LineWidth', 2);
plot(tsim, ynonlineaire(:,19), 'g--', 'LineWidth', 2);
plot(tsim, yLineaire.Data(:,1), 'b');
plot(tsim, yLineaire.Data(:,2), 'r');
plot(tsim, yLineaire.Data(:,3), 'g');
ylim([0 0.025])
legend('ZD (NL)', 'ZE (NL)', 'ZF (NL)', 'ZD (L)', 'ZE (L)', 'ZF (L)');grid on; grid minor;
subplot(212); 
hold on;title('Différence des distances ZA ZB ZC ZD ZE ZF');
plot(tsim, yLineaire.Data(:,1) - ynonlineaire(:,17), 'b');
plot(tsim, yLineaire.Data(:,2) - ynonlineaire(:,18), 'r');
plot(tsim, yLineaire.Data(:,3) - ynonlineaire(:,19), 'g');
legend('\Delta ZD', '\Delta ZE', '\Delta ZF');grid on; grid minor;
ylim([-0.01 0.01])

