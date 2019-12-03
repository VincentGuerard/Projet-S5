%%
%Identification des capteurs
%Aymen Mousli
%Projet
%18-11-2019 11:41

clear all
close all
clc

load('../DonneesExperimentales/capteur.mat')

%% Modele 1 - polynomiale + cos + exponetielle
phi1 = -voltage.^3;
phi2 = -voltage.^2;
phi3 = cos(voltage);
phi4 = -1./(exp(1./voltage));

Phi1Phi1 = sum(phi1.*phi1);
Phi2Phi2 = sum(phi2.*phi2);
Phi3Phi3 = sum(phi3.*phi3);
Phi4Phi4 = sum(phi4.*phi4);

Phi1Phi2 = sum(phi1.*phi2);
Phi1Phi3 = sum(phi1.*phi3);
Phi1Phi4 = sum(phi1.*phi4);
Phi2Phi3 = sum(phi2.*phi3);
Phi2Phi4 = sum(phi2.*phi4);
Phi3Phi4 = sum(phi3.*phi4);

A = [Phi1Phi1 Phi1Phi2 Phi1Phi3 Phi1Phi4;Phi1Phi2 Phi2Phi2 Phi2Phi3 Phi2Phi4;Phi1Phi3 Phi2Phi3 Phi3Phi3 Phi3Phi4;Phi1Phi4 Phi2Phi4 Phi3Phi4 Phi4Phi4];

Phi1Y = sum(phi1.*distance);
Phi2Y = sum(phi2.*distance);
Phi3Y = sum(phi3.*distance);
Phi4Y = sum(phi4.*distance);

YY= [Phi1Y;Phi2Y;Phi3Y;Phi4Y];

XX=inv(A)*YY;

yfin = XX(1)*phi1 + XX(2)*phi2 + XX(3)*phi3 + XX(4)*phi4;

figure
plot(voltage,distance,'.', 'MarkerSize', 8)
hold on
plot(voltage,yfin, 'Linewidth', 2)
title('Modèle 1 - polynomial + cos + exponentiel');
xlabel('Tension (V)');
ylabel('Distance (m)');
legend('Mesures', 'Fonction approximée');
grid on; grid minor;

%% Modele 2 - Sinus et cosinus
phi1 = cos(voltage);
phi2 = sin(2*voltage);
phi3 = cos(4*voltage);
phi4 = sin(6*voltage);
phi5 = cos(voltage);
phi6 = sin(voltage);

Phi1Phi1 = sum(phi1.*phi1);
Phi2Phi2 = sum(phi2.*phi2);
Phi3Phi3 = sum(phi3.*phi3);
Phi4Phi4 = sum(phi4.*phi4);
Phi5Phi5 = sum(phi5.*phi5);
Phi6Phi6 = sum(phi6.*phi6);

Phi1Phi2 = sum(phi1.*phi2);
Phi1Phi3 = sum(phi1.*phi3);
Phi1Phi4 = sum(phi1.*phi4);
Phi1Phi5 = sum(phi1.*phi5);
Phi1Phi6 = sum(phi1.*phi6);
Phi2Phi3 = sum(phi2.*phi3);
Phi2Phi5 = sum(phi2.*phi5);
Phi2Phi6 = sum(phi2.*phi6);
Phi2Phi4 = sum(phi2.*phi4);
Phi3Phi4 = sum(phi3.*phi4);
Phi3Phi5 = sum(phi3.*phi4);
Phi3Phi6 = sum(phi3.*phi4);
Phi4Phi5 = sum(phi4.*phi5);
Phi4Phi6 = sum(phi4.*phi6);
Phi5Phi6 = sum(phi5.*phi6);

A = [Phi1Phi1 Phi1Phi2 Phi1Phi3 Phi1Phi4 Phi1Phi5 Phi1Phi6;
    Phi1Phi2 Phi2Phi2 Phi2Phi3 Phi2Phi4 Phi2Phi5 Phi2Phi6;
    Phi1Phi3 Phi2Phi3 Phi3Phi3 Phi3Phi4 Phi3Phi5 Phi3Phi6;
    Phi1Phi4 Phi2Phi4 Phi3Phi4 Phi4Phi4 Phi4Phi5 Phi4Phi6;
    Phi1Phi5 Phi2Phi5 Phi3Phi5 Phi4Phi5 Phi5Phi5 Phi5Phi6;
    Phi1Phi6 Phi2Phi6 Phi3Phi6 Phi4Phi6 Phi5Phi6 Phi6Phi6];

Phi1Y = sum(phi1.*distance);
Phi2Y = sum(phi2.*distance);
Phi3Y = sum(phi3.*distance);
Phi4Y = sum(phi4.*distance);
Phi5Y = sum(phi5.*distance);
Phi6Y = sum(phi6.*distance);

YY= [Phi1Y;Phi2Y;Phi3Y;Phi4Y;;Phi5Y;Phi6Y];

XX=inv(A)*YY;

yfin = XX(1)*phi1 + XX(2)*phi2 + XX(3)*phi3 + XX(4)*phi4 + XX(5)*phi5 + XX(6)*phi6;

figure
plot(voltage,distance,'.', 'MarkerSize' ,8)
hold on
plot(voltage,yfin, 'Linewidth', 2)
title('Modèle 2 - Sinus et Cosinus');
xlabel('Tension (V)');
ylabel('Distance (m)');
legend('Mesures', 'Fonction approximée');
grid on; grid minor;

%% Modele 3 - Exponetielles et logarithmes
phi1 = exp(voltage);
phi2 = log(2*voltage);
phi3 = exp(4*voltage);
phi4 = log(6*voltage);
phi5 = exp(5*voltage);
phi6 = log(voltage);

Phi1Phi1 = sum(phi1.*phi1);
Phi2Phi2 = sum(phi2.*phi2);
Phi3Phi3 = sum(phi3.*phi3);
Phi4Phi4 = sum(phi4.*phi4);
Phi5Phi5 = sum(phi5.*phi5);
Phi6Phi6 = sum(phi6.*phi6);

Phi1Phi2 = sum(phi1.*phi2);
Phi1Phi3 = sum(phi1.*phi3);
Phi1Phi4 = sum(phi1.*phi4);
Phi1Phi5 = sum(phi1.*phi5);
Phi1Phi6 = sum(phi1.*phi6);
Phi2Phi3 = sum(phi2.*phi3);
Phi2Phi5 = sum(phi2.*phi5);
Phi2Phi6 = sum(phi2.*phi6);
Phi2Phi4 = sum(phi2.*phi4);
Phi3Phi4 = sum(phi3.*phi4);
Phi3Phi5 = sum(phi3.*phi4);
Phi3Phi6 = sum(phi3.*phi4);
Phi4Phi5 = sum(phi4.*phi5);
Phi4Phi6 = sum(phi4.*phi6);
Phi5Phi6 = sum(phi5.*phi6);

A = [Phi1Phi1 Phi1Phi2 Phi1Phi3 Phi1Phi4 Phi1Phi5 Phi1Phi6;
    Phi1Phi2 Phi2Phi2 Phi2Phi3 Phi2Phi4 Phi2Phi5 Phi2Phi6;
    Phi1Phi3 Phi2Phi3 Phi3Phi3 Phi3Phi4 Phi3Phi5 Phi3Phi6;
    Phi1Phi4 Phi2Phi4 Phi3Phi4 Phi4Phi4 Phi4Phi5 Phi4Phi6;
    Phi1Phi5 Phi2Phi5 Phi3Phi5 Phi4Phi5 Phi5Phi5 Phi5Phi6;
    Phi1Phi6 Phi2Phi6 Phi3Phi6 Phi4Phi6 Phi5Phi6 Phi6Phi6];

Phi1Y = sum(phi1.*distance);
Phi2Y = sum(phi2.*distance);
Phi3Y = sum(phi3.*distance);
Phi4Y = sum(phi4.*distance);
Phi5Y = sum(phi5.*distance);
Phi6Y = sum(phi6.*distance);

YY= [Phi1Y;Phi2Y;Phi3Y;Phi4Y;;Phi5Y;Phi6Y];

XX=inv(A)*YY;

yfin = XX(1)*phi1 + XX(2)*phi2 + XX(3)*phi3 + XX(4)*phi4 + XX(5)*phi5 + XX(6)*phi6;

figure
plot(voltage,distance, '.', 'MarkerSize', 8)
hold on
plot(voltage,yfin, 'Linewidth', 2)
title('Modèle 3 - Exponentielle et logarithme');
xlabel('Tension (V)');
ylabel('Distance (m)');
legend('Mesures', 'Fonction approximée');
grid on; grid minor;

%% Modele 4 - Sinh et cosh
phi1 = cosh(voltage);
phi2 = sinh(voltage);
phi3 = cosh(2*voltage);
phi4 = sinh(3*voltage);

Phi1Phi1 = sum(phi1.*phi1);
Phi2Phi2 = sum(phi2.*phi2);
Phi3Phi3 = sum(phi3.*phi3);
Phi4Phi4 = sum(phi4.*phi4);

Phi1Phi2 = sum(phi1.*phi2);
Phi1Phi3 = sum(phi1.*phi3);
Phi1Phi4 = sum(phi1.*phi4);
Phi2Phi3 = sum(phi2.*phi3);
Phi2Phi4 = sum(phi2.*phi4);
Phi3Phi4 = sum(phi3.*phi4);

A = [Phi1Phi1 Phi1Phi2 Phi1Phi3 Phi1Phi4;Phi1Phi2 Phi2Phi2 Phi2Phi3 Phi2Phi4;Phi1Phi3 Phi2Phi3 Phi3Phi3 Phi3Phi4;Phi1Phi4 Phi2Phi4 Phi3Phi4 Phi4Phi4];

Phi1Y = sum(phi1.*distance);
Phi2Y = sum(phi2.*distance);
Phi3Y = sum(phi3.*distance);
Phi4Y = sum(phi4.*distance);

YY= [Phi1Y;Phi2Y;Phi3Y;Phi4Y];

XX=inv(A)*YY;

yfin = XX(1)*phi1 + XX(2)*phi2 + XX(3)*phi3 + XX(4)*phi4;

figure
plot(voltage,distance,'.', 'MarkerSize', 8)
hold on
plot(voltage,yfin, 'LineWidth', 2)
title('Modèle 4 - Sinus hyperbolique et Cosinus hyperbolique');
xlabel('Tension (V)');
ylabel('Distance (m)');
legend('Mesures', 'Fonction approximée');
grid on; grid minor;
