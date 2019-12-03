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
% Ts = 5;
% t_des     = [0:1:10]'*Ts;
% x = [0.0106    0.0212    0.0317    0.0421    0.0522    0.0618    0.0687    0.0659    0.0580    0.0491    0.0398]';
% y = [0.0434    0.0646    0.0688    0.0610    0.0465    0.0303    0.0193    0.0235    0.0367    0.0514    0.0634]';
% x_des     = [t_des, x];
% y_des     = [t_des, y];
% z_des     = [t_des, [1 1 1 1  1  1 1 1 1 1 1]'*.015];
% tfin = 50;

% %Test trajectoire simple
% t_des     = [0:1:2]'*5;
% x_des     = [t_des, [0 0 0]'*0.05];
% y_des     = [t_des, [0 1 1]'*0.05];
% z_des     = [t_des, [1 1 1]'*.015];
% tfin = 15;

%Trejectoir référence
positionAB = [0 0.2 0.4 0.6 0.8 1;
              0 .3 .5 0.3 0.4 1]'*0.05;    
v = 0.01;
Ts = 0.2;
[Pi, Ltr, E, Vr, Traj, tt] = Trajectoire_func(positionAB,v,Ts);
x = Traj(:,1);
y = Traj(:,2);
t_des = [0:1:length(x)-1]'.*Ts;
x_des = [t_des, x];
y_des = [t_des, y];
z_des  = [t_des, ones(length(x),1)*.015];
tfin = tt+0.5+0.2;

% %Test hauteur
% t_des     = [0:1:2]'*5;
% x_des     = [t_des, [0 0 0]'*0.05];
% y_des     = [t_des, [0 0 0]'*0.05];
% z_des     = [t_des, [1 0.5 1]'*.015];
% tfin = 15;

%Test à une serie de mouvements
% t_des     = [0:1:9]'*5;

% x_des     = [t_des, [0 1 2 2 0 0 0 -1 1 0]'*0.05];
% y_des     = [t_des, [0 2 2 -2 0 0 0 2 1 1]'*0.05];
% z_des     = [t_des, [1 1 1 1 1 1 1 1 1 1]'*.015];
% tfin = 50;

%Tester Trajectoire 
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
numCompAngles = 1.0e09*[0.000230814504577 0.026170635469275 0.897433010674818 9.699372248008062];
denCompAngles = 1.0e07*[0.000000100000000 0.000823402569667 1.694979479335234 0];
numCompz = 1.0e11*[0.000031609546854 0.004304100455063 0.142047621936623 1.354942422178642];
denCompz = 1.0e05*[0.000010000000000 0.017380342819025 5.055406372627656 0];

%simulation
open_system('DYNctl_ver4_etud_obfusc')
set_param('DYNctl_ver4_etud_obfusc','AlgebraicLoopSolver','TrustRegion')
sim('DYNctl_ver4_etud_obfusc')

% Verification des violations
DetectionViolations

%% affichage

%Figure test trajectoire référence
figure;subplot(211);hold on
plot(ySystemeNonLineaire(:,2), ySystemeNonLineaire(:,1),'b');
plot(positionAB(:,1), positionAB(:,2), 'o')
legend('Px', 'Py')
hold off;
subplot(212); plot(tsim, ySystemeNonLineaire(:,4), 'b'); plot(tsim, ySystemeNonLineaire(:,3), 'r');
legend('v_x', 'v_y')
errx = positionAB(end,1)-ySystemeNonLineaire(end,2)
erry = positionAB(end,2)-ySystemeNonLineaire(end,1)

% %Figure test simple
% figure; hold on; plot((0:1:14), [0 0 0 0 0 1 1 1 1 1 0 0 0 0 0]*0.05, 'b', 'Linewidth', 2); plot((0:1:14), [0 0 0 0 0 0 0 0 0 0 1 1 1 1 1]*0.05, 'r', 'Linewidth', 2);
% plot(tsim, ySystemeNonLineaire(:,1),'b'); plot(tsim, ySystemeNonLineaire(:,2),'r');
% legend('xdes', 'ydes', 'Px', 'Py')
% hold off

%
figure; hold on; plot(tsim, Ax_des1, 'b', 'Linewidth', 2); plot(tsim, Ay_des1, 'r', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,8), 'b');
plot(tsim, ySystemeNonLineaire(:,9), 'r');
legend('Ax_d_e_s', 'Ay_d_e_s', 'Ax', 'Ay');grid on; grid minor;
hold off

figure; hold on; plot((0:1:14), [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]*0.015, 'b--', 'Linewidth', 2); plot(tsim, ySystemeNonLineaire(:,10),'b')
legend('z_d_e_s', 'Pz')
hold off

% Comparaison des tension Va, Vb, Vc
figure
subplot(211); 
plot(tsim, ynonlineaire(:,23), 'b--');hold on;title('Tension Va Vb Vc');
plot(tsim, ynonlineaire(:,24), 'r--');
plot(tsim, ynonlineaire(:,25), 'g--');
plot(tsim, ySystemeNonLineaire(:,17), 'b');
plot(tsim, ySystemeNonLineaire(:,18), 'r');
plot(tsim, ySystemeNonLineaire(:,19), 'g');
legend('Va prof', 'Vb prof', 'Vc prof', 'Va', 'Vb', 'Vc');grid on; grid minor;
subplot(212);
plot(tsim, ynonlineaire(:,23) - ySystemeNonLineaire(:,17), 'b'); hold on; title('Différence des tension Va Vb Vc');
plot(tsim, ynonlineaire(:,24) - ySystemeNonLineaire(:,18), 'r');
plot(tsim, ynonlineaire(:,25) - ySystemeNonLineaire(:,19), 'g');
legend('Diff Va', 'Diff Vb', 'Diff Vc');grid on; grid minor;

% % Comparaison des courants Ia, Ib, Ic
% figure
% subplot(211); 
% plot(tsim, ynonlineaire(:,11), '--');hold on;title('Courant Ia Ib Ic des deux modèles');
% plot(tsim, ynonlineaire(:,12), '--');
% plot(tsim, ynonlineaire(:,13), '--');
% plot(tsim, ySystemeNonLineaire(:,5));
% plot(tsim, ySystemeNonLineaire(:,6));
% plot(tsim, ySystemeNonLineaire(:,7));
% legend('Ia prof', 'Ib prof', 'Ic prof', 'Ia', 'Ib', 'Ic');grid on; grid minor;
% subplot(212); 
% plot(tsim, ySystemeNonLineaire(:,5) - ynonlineaire(:,11)); hold on;title('Différence des courants Ia Ib Ic');
% plot(tsim, ySystemeNonLineaire(:,6) - ynonlineaire(:,12));
% plot(tsim, ySystemeNonLineaire(:,7) - ynonlineaire(:,13));
% legend('Diff Ia', 'Diff Ib', 'Diff Ic');grid on; grid minor;

% % Comparaison des Forces Fa, Fb, Fc
% figure
% subplot(211); 
% plot(tsim, ynonlineaire(:,20), '--');hold on;title('Forces Fa Fb Fc des deux modèles');
% plot(tsim, ynonlineaire(:,21), '--');
% plot(tsim, ynonlineaire(:,22), '--');
% plot(tsim, ySystemeNonLineaire(:,14));
% plot(tsim, ySystemeNonLineaire(:,15));
% plot(tsim, ySystemeNonLineaire(:,16));
% legend('Fa prof', 'Fb prof', 'Fc prof', 'Fa', 'Fb', 'Fc');grid on; grid minor;
% ylim([-5 5]);
% subplot(212); 
% plot(tsim, ySystemeNonLineaire(:,14) - ynonlineaire(:,20)); hold on;title('Différence des forces Fa Fb Fc');
% plot(tsim, ySystemeNonLineaire(:,15) - ynonlineaire(:,21));
% plot(tsim, ySystemeNonLineaire(:,16) - ynonlineaire(:,22));
% legend('Diff Fa', 'Diff Fb', 'Diff Fc');grid on; grid minor;

% Comparaison des angles de la plaque phi et theta
figure
plot(tsim, ynonlineaire(:,1), 'b--', 'Linewidth', 2);hold on;title('Angles de la plaque');
plot(tsim, ynonlineaire(:,2), 'r--', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,9), 'b');
plot(tsim, ySystemeNonLineaire(:,8), 'r');
legend('Ax prof', 'Ay prof', 'Ax', 'Ay');grid on; grid minor;

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
subplot(212); 
plot(tsim, ySystemeNonLineaire(:,2) - ynonlineaire(:,7), 'b', 'Linewidth', 2); hold on;title('Différence de la position Px Py et la vitesse Vx Vy');
plot(tsim, ySystemeNonLineaire(:,1) - ynonlineaire(:,8), 'r', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,4) - ynonlineaire(:,9), 'g', 'Linewidth', 2);
plot(tsim, ySystemeNonLineaire(:,3) - ynonlineaire(:,10), 'y', 'Linewidth', 2);
legend('Diff Px', 'Diff Py', 'Diff Vx', 'Diff Vy');grid on; grid minor;

% %% Comparaison des distances
% figure
% subplot(211); 
% plot(tsim, ynonlineaire(:,14), '--');
% hold on;title('Distance zA zB zC zD zE zF des deux modèles');
% plot(tsim, ynonlineaire(:,15), '--');
% plot(tsim, ynonlineaire(:,16), '--');
% plot(tsim, ynonlineaire(:,17), '--');
% plot(tsim, ynonlineaire(:,18), '--');
% plot(tsim, ynonlineaire(:,19), '--');
% plot(tsim, ySystemeNonLineaire(:,20));
% plot(tsim, ySystemeNonLineaire(:,21));
% plot(tsim, ySystemeNonLineaire(:,22));
% plot(tsim, ySystemeNonLineaire(:,23));
% plot(tsim, ySystemeNonLineaire(:,24));
% plot(tsim, ySystemeNonLineaire(:,25));
% legend('zA prof', 'zB prof', 'zC prof', 'zD prof', 'zE prof', 'zF prof', 'zA', 'zB', 'zC', 'zD', 'zE', 'zF');grid on; grid minor;
% subplot(212); 
% plot(tsim, ySystemeNonLineaire(:,20) - ynonlineaire(:,14)); hold on;title('Différence des distances zA zB zC zD zE zF');
% plot(tsim, ySystemeNonLineaire(:,21) - ynonlineaire(:,15));
% plot(tsim, ySystemeNonLineaire(:,22) - ynonlineaire(:,16));
% plot(tsim, ySystemeNonLineaire(:,23) - ynonlineaire(:,17));
% plot(tsim, ySystemeNonLineaire(:,24) - ynonlineaire(:,18));
% plot(tsim, ySystemeNonLineaire(:,25) - ynonlineaire(:,19));
% legend('Diff zA', 'Diff zB', 'Diff zC', 'Diff zD', 'Diff zE', 'Diff zF');grid on; grid minor;


