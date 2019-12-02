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
%bancessai_ini  %faites tous vos calculs de modele ici
%run('../../Identification/IdentificationActionneur.m');close all;
run('../../Identification/Projet_s5.m');close all;
bancEssaiConstantes

%Calcul des compensateurs
%iniCTL_ver4    %Calculez vos compensateurs ici
Tdef = [YD -XD 1;
        YE -XE 1;
        YF -XF 1];
TABC = [YA YB YC;
        -XA -XB -XC;
        1 1 1];
%simulation
open_system('DYNctl_ver4_etud_obfusc')
set_param('DYNctl_ver4_etud_obfusc','AlgebraicLoopSolver','LineSearch')
sim('DYNctl_ver4_etud_obfusc')

%% affichage

% Comparaison des tension Va, Vb, Vc
figure
subplot(211); 
plot(tsim, ynonlineaire(:,23), '--');hold on;title('Tension Va Vb Vc des deux modèles');
plot(tsim, ynonlineaire(:,24), '--');
plot(tsim, ynonlineaire(:,25), '--');
plot(tsim, ySystemeNonLineaire(:,17));
plot(tsim, ySystemeNonLineaire(:,18));
plot(tsim, ySystemeNonLineaire(:,19));
legend('Va prof', 'Vb prof', 'Vc prof', 'Va', 'Vb', 'Vc');grid on; grid minor;
subplot(212);
plot(tsim, ynonlineaire(:,23) - ySystemeNonLineaire(:,17)); hold on; title('Différence des tension Va Vb Vc');
plot(tsim, ynonlineaire(:,24) - ySystemeNonLineaire(:,18));
plot(tsim, ynonlineaire(:,25) - ySystemeNonLineaire(:,19));
legend('Diff Va', 'Diff Vb', 'Diff Vc');grid on; grid minor;

% Comparaison des courants Ia, Ib, Ic
figure
subplot(211); 
plot(tsim, ynonlineaire(:,11), '--');hold on;title('Courant Ia Ib Ic des deux modèles');
plot(tsim, ynonlineaire(:,12), '--');
plot(tsim, ynonlineaire(:,13), '--');
plot(tsim, ySystemeNonLineaire(:,5));
plot(tsim, ySystemeNonLineaire(:,6));
plot(tsim, ySystemeNonLineaire(:,7));
legend('Ia prof', 'Ib prof', 'Ic prof', 'Ia', 'Ib', 'Ic');grid on; grid minor;
subplot(212); 
plot(tsim, ySystemeNonLineaire(:,5) - ynonlineaire(:,11)); hold on;title('Différence des courants Ia Ib Ic');
plot(tsim, ySystemeNonLineaire(:,6) - ynonlineaire(:,12));
plot(tsim, ySystemeNonLineaire(:,7) - ynonlineaire(:,13));
legend('Diff Ia', 'Diff Ib', 'Diff Ic');grid on; grid minor;

% Comparaison des Forces Fa, Fb, Fc
figure
subplot(211); 
plot(tsim, ynonlineaire(:,20), '--');hold on;title('Forces Fa Fb Fc des deux modèles');
plot(tsim, ynonlineaire(:,21), '--');
plot(tsim, ynonlineaire(:,22), '--');
plot(tsim, ySystemeNonLineaire(:,14));
plot(tsim, ySystemeNonLineaire(:,15));
plot(tsim, ySystemeNonLineaire(:,16));
legend('Fa prof', 'Fb prof', 'Fc prof', 'Fa', 'Fb', 'Fc');grid on; grid minor;
ylim([-5 5]);
subplot(212); 
plot(tsim, ySystemeNonLineaire(:,14) - ynonlineaire(:,20)); hold on;title('Différence des forces Fa Fb Fc');
plot(tsim, ySystemeNonLineaire(:,15) - ynonlineaire(:,21));
plot(tsim, ySystemeNonLineaire(:,16) - ynonlineaire(:,22));
legend('Diff Fa', 'Diff Fb', 'Diff Fc');grid on; grid minor;

% Comparaison des angles de la plaque phi et theta
figure
plot(tsim, ynonlineaire(:,1), '--');hold on;title('Angle des deux modèles');
plot(tsim, ynonlineaire(:,2), '--');
plot(tsim, ynonlineaire(:,4), '--');
plot(tsim, ynonlineaire(:,5), '--');
plot(tsim, ySystemeNonLineaire(:,9));
plot(tsim, ySystemeNonLineaire(:,8));
plot(tsim, ySystemeNonLineaire(:,12));
plot(tsim, ySystemeNonLineaire(:,13));
legend('Ax prof', 'Ay prof', 'Wx prof', 'Wy prof', 'Ax', 'Ay', 'Wx', 'Wy');grid on; grid minor;

% Comparaison de la position de la bille sur la plaque Px et Py
figure
subplot(211); 
plot(tsim, ynonlineaire(:,7), 'b--', 'Linewidth', 2);hold on;title('Position Px Py Vx Vy des deux modèles');
plot(tsim, ynonlineaire(:,8), 'r--', 'Linewidth', 2);
plot(tsim, ynonlineaire(:,9), 'g--', 'Linewidth', 2);
plot(tsim, ynonlineaire(:,10), 'y--', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,2),'b');
plot(tsim, ySystemeNonLineaire(:,1),'r');
plot(tsim, ySystemeNonLineaire(:,4), 'g');
plot(tsim, ySystemeNonLineaire(:,3), 'y');
legend('Px prof', 'Py prof', 'Vx prof', 'Vy prof', 'Px', 'Py', 'Vx', 'Vy');grid on; grid minor;
ylim([-0.06 0.06]);
subplot(212); 
plot(tsim, ySystemeNonLineaire(:,2) - ynonlineaire(:,7), 'b', 'Linewidth', 2); hold on;title('Différence de la position Px Py et la vitesse Vx Vy');
plot(tsim, ySystemeNonLineaire(:,1) - ynonlineaire(:,8), 'r', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,4) - ynonlineaire(:,9), 'g', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,3) - ynonlineaire(:,10), 'y', 'Linewidth', 2);
legend('Diff Px', 'Diff Py', 'Diff Vx', 'Diff Vy');grid on; grid minor;
ylim([-0.06 0.06])

%% Comparaison des distances
figure
subplot(211); 
plot(tsim, ynonlineaire(:,14), '--');
hold on;title('Distance zA zB zC zD zE zF des deux modèles');
plot(tsim, ynonlineaire(:,15), '--');
plot(tsim, ynonlineaire(:,16), '--');
plot(tsim, ynonlineaire(:,17), '--');
plot(tsim, ynonlineaire(:,18), '--');
plot(tsim, ynonlineaire(:,19), '--');
plot(tsim, ySystemeNonLineaire(:,20));
plot(tsim, ySystemeNonLineaire(:,21));
plot(tsim, ySystemeNonLineaire(:,22));
plot(tsim, ySystemeNonLineaire(:,23));
plot(tsim, ySystemeNonLineaire(:,24));
plot(tsim, ySystemeNonLineaire(:,25));
legend('zA prof', 'zB prof', 'zC prof', 'zD prof', 'zE prof', 'zF prof', 'zA', 'zB', 'zC', 'zD', 'zE', 'zF');grid on; grid minor;
subplot(212); 
plot(tsim, ySystemeNonLineaire(:,20) - ynonlineaire(:,14)); hold on;title('Différence des distances zA zB zC zD zE zF');
plot(tsim, ySystemeNonLineaire(:,21) - ynonlineaire(:,15));
plot(tsim, ySystemeNonLineaire(:,22) - ynonlineaire(:,16));
plot(tsim, ySystemeNonLineaire(:,23) - ynonlineaire(:,17));
plot(tsim, ySystemeNonLineaire(:,24) - ynonlineaire(:,18));
plot(tsim, ySystemeNonLineaire(:,25) - ynonlineaire(:,19));
legend('Diff zA', 'Diff zB', 'Diff zC', 'Diff zD', 'Diff zE', 'Diff zF');grid on; grid minor;


