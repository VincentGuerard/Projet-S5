%% Identification de la formule de l'actionneur
clear all
close all
clc

load("..\DonneesExperimentales\Fe_attraction.mat")
load("..\DonneesExperimentales\Fs.mat")

%% Identification de la fonction Fs avec la méthode des Moindres Carrées
precision = 1e-3;
xn = z_pos;
yn = Fs;
parametres = [0.1
              0.1;
              0.1;
              0.1];

for i = 0:1:50
    Jr = [-1./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          -xn./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (-xn.^2)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (-xn.^3)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2];

    Jrinv = pinv(Jr);
    r = yn - (-1./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)));
    r2 = r.^2;
    
    parametres = parametres - (Jrinv * r);
    
    if abs(sum(r2)) < precision
        break
    end
end
parametresFs = parametres;

gx = -1./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1));
x = 0:0.0001:0.035;
gxc = -1./(parametres(4).*x.^3 + parametres(3).*x.^2 + parametres(2).*x + parametres(1));
N = length(xn);
RMS = sqrt(1/N*sum((gx-yn).^2));
yprime = 1/N*sum(yn);
SSres = sum((yn - gx).^2);
SST = sum((yn - yprime).^2);
R = 1 - SSres/SST;

figure
plot(x,gxc); grid on; grid minor; title('Fs en fonction de z_{pos}'); xlabel('Position (m)'); ylabel('Force (N)');
hold on
scatter(xn,yn);
disp(['Pour la fonction Fs, la valeur de R^2 est de: ' num2str(R) ' et la valeur RMS est de: ' num2str(RMS)]);
fprintf('\n');

%% Identification de la fonction Fe avec la méthode des Moindres Carrées
precision = 1e-1;
be1 = 13.029359254409743;

% 1 A
xn = z_m1A;
yn = Fe_m1A;
ik = -1;
parametres = [0.1
              0.1;
              0.1;
              0.1];

for i = 0:1:50
    Jr = [((ik^2+be1*abs(ik))*sign(ik))./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (((ik^2+be1*abs(ik))*sign(ik)).*xn)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (((ik^2+be1*abs(ik))*sign(ik)).*xn.^2)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (((ik^2+be1*abs(ik))*sign(ik)).*xn.^3)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2];

    Jrinv = pinv(Jr);
    r = yn - (((ik^2+be1*abs(ik))*sign(ik)) ./ (parametres(1) + parametres(2).*xn + parametres(3).*xn.^2 + parametres(4).*xn.^3));
    r2 = r.^2;
    
    parametres = parametres - (Jrinv * r);
    
    if abs(sum(r2)) < precision
        break
    end
end
parametres1A = parametres;

gx = ((ik^2+be1*abs(ik))*sign(ik)) ./ (parametres(1) + parametres(2).*xn + parametres(3).*xn.^2 + parametres(4).*xn.^3);
x = 0:0.0001:0.035;
gxc = ((ik^2+be1*abs(ik))*sign(ik)) ./ (parametres(1) + parametres(2).*x + parametres(3).*x.^2 + parametres(4).*x.^3);
N = length(xn);
RMS = sqrt(1/N*sum((gx-yn).^2));
yprime = 1/N*sum(yn);
SSres = sum((yn - gx).^2);
SST = sum((yn - yprime).^2);
R = 1 - SSres/SST;

figure
plot(x,gxc); grid on; grid minor; title('Fe avec courant de 1A en fonction de z_{pos}'); xlabel('Position (m)'); ylabel('Fe (N)');
hold on
scatter(xn,yn);
disp(['Pour la fonction Fe avec courant de 1A, la valeur de R^2 est de: ' num2str(R) ' et la valeur RMS est de: ' num2str(RMS)]);
fprintf('\n');


% 2 A
xn = z_m2A;
yn = Fe_m2A;
ik = -2;
parametres = [0.1
              0.1;
              0.1;
              0.1];

for i = 0:1:1000
    Jr = [((ik^2+be1*abs(ik))*sign(ik))./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (((ik^2+be1*abs(ik))*sign(ik)).*xn)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (((ik^2+be1*abs(ik))*sign(ik)).*xn.^2)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2 ...
          (((ik^2+be1*abs(ik))*sign(ik)).*xn.^3)./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1)).^2];

    Jrinv = pinv(Jr);
    r = yn - (((ik^2+be1*abs(ik))*sign(ik)) ./ (parametres(1) + parametres(2).*xn + parametres(3).*xn.^2 + parametres(4).*xn.^3));
    r2 = r.^2;
    
    parametres = parametres - (Jrinv * r);
    
    if abs(sum(r2)) < precision
        break
    end
end
parametres2A = parametres;

gx = ((ik^2+be1*abs(ik))*sign(ik)) ./ (parametres(1) + parametres(2).*xn + parametres(3).*xn.^2 + parametres(4).*xn.^3);
x = 0:0.0001:0.035;
gxc = ((ik^2+be1*abs(ik))*sign(ik)) ./ (parametres(1) + parametres(2).*x + parametres(3).*x.^2 + parametres(4).*x.^3);
N = length(xn);
RMS = sqrt(1/N*sum((gx-yn).^2));
yprime = 1/N*sum(yn);
SSres = sum((yn - gx).^2);
SST = sum((yn - yprime).^2);
R = 1 - SSres/SST;

figure
plot(x,gxc); grid on; grid minor; title('Fe avec courant de 2A en fonction de z_{pos}'); xlabel('Position (m)'); ylabel('Fe (N)');
hold on
scatter(xn,yn);
disp(['Pour la fonction Fe avec courant de 2A, la valeur de R^2 est de: ' num2str(R) ' et la valeur RMS est de: ' num2str(RMS)]);
fprintf('\n');

% Mapping vers les parametres du simulink
as0 = parametresFs(1);
as1 = parametresFs(2);
as2 = parametresFs(3);
as3 = parametresFs(4);

ae0 = parametres1A(1);
ae1 = parametres1A(2);
ae2 = parametres1A(3);
ae3 = parametres1A(4);


