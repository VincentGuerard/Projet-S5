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
run('../Projet_S5_decouple.m');

%Calcul des compensateurs
%iniCTL_ver4    %Calculez vos compensateurs ici
Compensateur_angles_v1

%simulation
open_system('DYNctl_ver4_etud_obfusc')
set_param('DYNctl_ver4_etud_obfusc','AlgebraicLoopSolver','LineSearch')
sim('DYNctl_ver4_etud_obfusc')

%% affichage

% Comparaison de la position de la bille sur la plaque Px et Py
figure
subplot(211); 
plot(tsim, ynonlineaire(:,7), 'b--', 'Linewidth', 2);hold on;title('Position Px Py Vx Vy nonlinéaire (pointillé) vs Linéaire (ligne pleine)');
plot(tsim, ynonlineaire(:,8), 'r--', 'Linewidth', 2);
plot(tsim, ynonlineaire(:,9), 'g--', 'Linewidth', 2);
plot(tsim, ynonlineaire(:,10), 'y--', 'Linewidth', 2);
plot(tsim, ySphereDecouple.Data(:,1), 'b');
plot(tsim, ySphereDecouple.Data(:,2), 'r');
plot(tsim, ySphereDecouple.Data(:,3), 'g');
plot(tsim, ySphereDecouple.Data(:,4), 'y');
legend('Px prof', 'Py prof', 'Vx prof', 'Vy prof', 'Px', 'Py', 'Vx', 'Vy');grid on; grid minor;
ylim([-0.06 0.06]);
subplot(212);
plot(tsim, ySphereDecouple.Data(:,1) - ynonlineaire(:,7), 'b', 'Linewidth', 2); hold on;title('Différence de la position Px Py et la vitesse Vx Vy');
plot(tsim, ySphereDecouple.Data(:,2) - ynonlineaire(:,8), 'r', 'Linewidth', 2);
plot(tsim, ySphereDecouple.Data(:,3) - ynonlineaire(:,9), 'g', 'Linewidth', 2);
plot(tsim, ySphereDecouple.Data(:,4) - ynonlineaire(:,10), 'y', 'Linewidth', 2);
legend('Diff Px', 'Diff Py', 'Diff Vx', 'Diff Vy');grid on; grid minor;

% Comparaison des distances
figure
subplot(211);
hold on;title('Distance zD zE zF des deux modèles');
plot(tsim, ynonlineaire(:,17), '--');
plot(tsim, ynonlineaire(:,18), '--');
plot(tsim, ynonlineaire(:,19), '--');
plot(tsim, yPlaqueDecouple.Data(:,1));
plot(tsim, yPlaqueDecouple.Data(:,2));
plot(tsim, yPlaqueDecouple.Data(:,3));
legend('zD prof', 'zE prof', 'zF prof', 'zD', 'zE', 'zF');grid on; grid minor;
subplot(212); 
hold on;title('Différence des distances zD zE zF');
plot(tsim, yPlaqueDecouple.Data(:,1) - ynonlineaire(:,17));
plot(tsim, yPlaqueDecouple.Data(:,2) - ynonlineaire(:,18));
plot(tsim, yPlaqueDecouple.Data(:,3) - ynonlineaire(:,19));
legend('Diff zD', 'Diff zE', 'Diff zF');grid on; grid minor;

% Comparaison des tenison
figure
subplot(211);
hold on;title('Tension Va Vb Vc des deux modèles');
plot(tsim, ynonlineaire(:,23), '--');
plot(tsim, ynonlineaire(:,24), '--');
plot(tsim, ynonlineaire(:,25), '--');
plot(tsim, Tension.Data(:,1));
plot(tsim, Tension.Data(:,2));
plot(tsim, Tension.Data(:,3));
legend('Va prof', 'Vb prof', 'Vc prof', 'Va', 'Vb', 'Vc');grid on; grid minor;
subplot(212);
hold on;title('Différence des distances Va Vb Vc');
plot(tsim, Tension.Data(:,1) - ynonlineaire(:,23));
plot(tsim, Tension.Data(:,2) - ynonlineaire(:,24));
plot(tsim, Tension.Data(:,3) - ynonlineaire(:,25));
legend('Diff Va', 'Diff Vb', 'Diff Vc');grid on; grid minor;

% Comparaison des angles
figure
subplot(211);
hold on;title('Angles et position des deux modèles');
plot(tsim, ynonlineaire(:,1), '--');
plot(tsim, ynonlineaire(:,2), '--');
plot(tsim, ynonlineaire(:,3), '--');
plot(tsim, Angles.Data(:,1));
plot(tsim, Angles.Data(:,2));
plot(tsim, Angles.Data(:,3));
legend('Phi prof', 'Theta prof', 'z prof', 'Phi', 'Theta', 'z');grid on; grid minor;
subplot(212); 
hold on;title('Différence de Phi Theta z');
plot(tsim, Angles.Data(:,1) - ynonlineaire(:,1));
plot(tsim, Angles.Data(:,2) - ynonlineaire(:,2));
plot(tsim, Angles.Data(:,3) - ynonlineaire(:,3));
legend('Diff Phi', 'Diff Theta', 'Diff z');grid on; grid minor;

