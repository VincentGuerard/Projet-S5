close all
clear all
clc

% Position � l'�quilibre de la sph�re (pour tests statiques)
sig = 0.0;         % Pr�sence (1) ou non (0) de la sph�re
xSeq = 0.000;      % Position x de la sph�re � l'�quilibre en metres
ySeq = 0.000;      % Position y de la sph�re � l'�quilibre en metres

%Point d'op�ration choisi pour la plaque
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
% run('../../Identification/IdentificationActionneur.m')
run('../../Identification/Projet_s5.m');close all;

%Calcul des compensateurs
%iniCTL_ver4    %Calculez vos compensateurs ici

%simulation
open_system('DYNctl_ver4_etud_obfusc')
set_param('DYNctl_ver4_etud_obfusc','AlgebraicLoopSolver','LineSearch')
sim('DYNctl_ver4_etud_obfusc')

%% affichage
% titre = ["Ax";"Ay";"Pz";"Wx";"Wy";"Vz";"Px";"Py";"Vx";"Vy";"Ia";"Ib";"Ic";"zA";"zB";"zC";"zD";"zE";"zF";"Fa";"Fb";"Fc";"Va";"Vb";"Vc"];
% for i = 1:1:25
%     figure
%     plot(tsim, ynonlineaire(:,i));title(titre(i));
% end
%trajectoires

% Comparaison des tension Va, Vb, Vc
figure
subplot(211)
plot(tsim, ynonlineaire(:,23));hold on;
plot(tsim, ynonlineaire(:,24));
plot(tsim, ynonlineaire(:,25));
plot(tsim, ySystemeNonLineaire(:,17));
plot(tsim, ySystemeNonLineaire(:,18));
plot(tsim, ySystemeNonLineaire(:,19));
subplot(212)
plot(tsim, ynonlineaire(:,23) - ySystemeNonLineaire(:,17)); hold on;
plot(tsim, ynonlineaire(:,24) - ySystemeNonLineaire(:,18));
plot(tsim, ynonlineaire(:,25) - ySystemeNonLineaire(:,19));

% Comparaison des courants Ia, Ib, Ic
figure
subplot(211)
plot(tsim, ynonlineaire(:,11));hold on;
plot(tsim, ynonlineaire(:,12));
plot(tsim, ynonlineaire(:,13));
plot(tsim, ySystemeNonLineaire(:,7));
plot(tsim, ySystemeNonLineaire(:,8));
plot(tsim, ySystemeNonLineaire(:,9));
subplot(212)
plot(tsim, ynonlineaire(:,11) - ySystemeNonLineaire(:,7)); hold on;
plot(tsim, ynonlineaire(:,12) - ySystemeNonLineaire(:,8));
plot(tsim, ynonlineaire(:,13) - ySystemeNonLineaire(:,9));

% Comparaison des Forces Fa, Fb, Fc

% Comparaison des angles de la plaque phi et theta
figure
plot(tsim, ynonlineaire(:,1));hold on;title('Angle des deux mod�les');
plot(tsim, ynonlineaire(:,2));
plot(tsim, ySystemeNonLineaire(:,9));
plot(tsim, ySystemeNonLineaire(:,8));
legend('Ax prof', 'Ay prof', 'Ax', 'Ay');grid on; grid minor;
ylim([-0.05 0.05]);

% Comparaison de la position de la bille sur la plaque Px et Py

