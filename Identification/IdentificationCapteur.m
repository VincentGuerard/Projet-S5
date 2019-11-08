%% Identification de la formule des capteurs
clear all
close all
clc

load("..\DonneesExperimentales\capteur.mat")

%% Modele 1 - Fourrier

x = voltage;
y = distance;

N = length(x);
w = (2*pi);

P1 = ones(length(x'),1);
P2 = cos(w*x);
P3 = sin(w*x);
P4 = cos(2*w*x);
P5 = sin(2*w*x);
P6 = cos(3*w*x);
P7 = sin(3*w*x);
P8 = cos(4*w*x);
P9 = sin(4*w*x);

P = [P1 P2 P3 P4 P5 P6 P7 P8 P9];
Ppinve = pinv(P);
A = Ppinve*y;

subplot(211)
scatter(x, y, 20);
hold on;
gx = A(1) + A(2)*cos(w*x) + A(3)*sin(w*x) + A(4)*cos(2*w*x) + A(5)*sin(2*w*x) + A(6)*cos(3*w*x) ...
     + A(7)*sin(3*w*x) + A(8)*cos(4*w*x) + A(9)*sin(4*w*x);
plot(x, gx);grid on; grid minor; title('Approximation par séries de Fourrier');xlabel('Distance (m)');ylabel('Tension (V)')

subplot(212)
difference = gx - y;
stem(x, difference);grid on; grid minor; title('Différence entre données et fonction');xlabel('Distance (m)');ylabel('Tension (V)');

RMS = sqrt(1/N*sum((gx-y).^2));
yprime = 1/N*sum(y);
R = sum((gx-yprime).^2)/sum((y-yprime).^2);

disp(['Pour l''approximation par série de Fourrier, l''erreur RMS est: ' num2str(RMS, '%.12f') ' et le facteur R est: ' num2str(R, '%.12f')]);
disp(['Avec ' num2str(length(A)) ' coefficients']);
fprintf('\n')

%% Méthode 2 - Forme exponentielle

x = voltage;
y = distance;

Y = log(y);

P1 = ones(length(x),1);
P2 = x;
P = [P1 P2];
Ppinve = pinv(P);
A = Ppinve*Y;
alpha1 = exp(A(1))
beta1 = A(2)

gx = beta1*exp(beta1.*x);

y = y - gx;
Y = log(y);
P1 = ones(length(x),1);
P2 = x;
P = [P1 P2];
Ppinve = pinv(P);
A1 = Ppinve*Y;
alpha2 = exp(A1(1))
beta2 = A1(2)

gx = alpha1*exp(beta1*x) - alpha2*exp(beta2*x);

figure
scatter(x,y, 20)
hold on
plot(x,gx)

RMS = sqrt(1/N*sum((gx-y).^2))















