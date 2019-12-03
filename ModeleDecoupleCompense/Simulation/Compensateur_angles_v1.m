%% Compensateur Angles
close all

% Specifications
Mp          = 0.05;
Ts          = 0.03;
Tp          = 0.025;
Tr          = 0.02;
erp_ech     = 0;
erp_rampe   = 1;

% Traduire les critères de performance
[p_des wn zeta] = criteresPerformanceTemporel(Mp, Ts, Tp, Tr, 0)

% Ajout d'un avance de phase pour réglé le régime transitoire
FTBO    = FTplaque(1,1);
Ga      = avPhase(p_des, FTBO, (15/360)*2*pi, 2);
FTBF    = FTBO * Ga;

% Affichage de l'atteinte des pôles désirés
figure;
rlocus(FTBF); hold on;
plot(p_des, 'p');
p = rlocus(FTBF, 1);
plot(real(p), imag(p), 's'); hold off
xlim([-1000 200]); ylim([-600 600])
title('Lieu des racines de Angles/V_{angles} avec AvPh')

%Vérification des performances
disp('Performances FTplaque des angles compensée avec AvPh:')
stepinfo(feedback(FTBF, 1))
t = 0:0.00001:0.25;
u = ones(length(t),1);

figure
lsim(feedback(FTBF,1),u,t)

% Ajout d'un PI pour réglé le regime permanent
[num den] = tfdata(FTBF, 'v');
z = real(p_des(1))/5;
numPI = [1 -z];
denPI = [1 0];
Kp = 1/abs((p_des(1) - z)/p_des(1) * (polyval(num, p_des(1))/polyval(den, p_des(1))));
Gpi = tf(Kp*[1 -z], [1 0]);

FTBF = FTBF * Gpi;

% Affichage de l'atteinte des pôles désirés
figure;
rlocus(FTBF); hold on;
plot(p_des, 'p');
p = rlocus(FTBF, 1);
plot(real(p), imag(p), 's'); hold off
xlim([-1000 200]); ylim([-600 600])
title('Lieu des racines de Angles / V_{angles} avec AvPh + PI')

%Vérification des performances
disp('Performances FTplaque des angles compensée avec AvPh + PI:')
stepinfo(feedback(FTBF, 1))
t = 0:0.00001:0.25;
u = ones(length(t),1);

figure
lsim(feedback(FTBF,1),u,t)

% Ajustement du gain
K = 1.5;
p = rlocus(FTBF, K);
angles = angle(p);
freq_n = real(p)./cos(angles);
while freq_n(2) < 1000 && freq_n(3) < 1000 && freq_n(4) < 1000 && freq_n(5) < 1000 && freq_n(6) < 1000 && K < 10
    K = K + 0.001;
    p = rlocus(FTBF, K);
    angles = angle(p);
    freq_n = real(p)./cos(angles);
end 
figure
rlocus(K * FTBF); hold on
p = rlocus(FTBF, K);
plot(real(p), imag(p), 's')

% Affichage de l'atteinte des pôles désirés
figure;
rlocus(K*FTBF); hold on;
plot(p_des, 'p');
p = rlocus(K*FTBF, 1);
plot(real(p), imag(p), 's'); hold off
xlim([-1000 200]); ylim([-600 600])
title('Lieu des racines de Angles / V_{angles} avec AvPh + PI')

%Vérification des performances
disp('Performances FTplaque des angles compensée avec AvPh + PI:')
stepinfo(feedback(K*FTBF, 1))
t = 0:0.00001:0.25;
u = ones(length(t),1);

figure
lsim(feedback(K*FTBF,1),u,t)

FTBF = FTBF * K
Comp = K * Ga * Gpi

% Vérification des critères de sécurités
figure 
margin(FTBF)

figure
nyquist(FTBF)
ylim([-3 3])











