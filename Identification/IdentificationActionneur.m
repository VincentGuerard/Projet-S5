%% Identification de la formule de l'actionneur
clear all
close all
clc

load("..\DonneesExperiemntales\Fe_attraction.mat")
load("..\DonneesExperiemntales\Fs.mat")

%% Identification de la fonction Fs avec la méthode de Gauss-Newton
precision = 1e-3;
xn = z_pos;
yn = Fs;

syms y b1 b2 x a1 a2 a3 a4
f = y - -1 / (a1 + a2*x + a3*x^2 + a4*x^3);
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

gx = -1./(parametres(4).*xn.^3 + parametres(3).*xn.^2 + parametres(2).*xn + parametres(1));
x = 0:0.0001:0.035;
gxc = -1./(parametres(4).*x.^3 + parametres(3).*x.^2 + parametres(2).*x + parametres(1));
N = length(xn);
RMS = sqrt(1/N*sum((gx-yn).^2));
yprime = 1/N*sum(yn);
R = sum((gx-yprime).^2)/sum((yn-yprime).^2);

figure
plot(x,gxc); grid on; grid minor; title('Fs en fonction de z_{pos}'); xlabel('Position (m)'); ylabel('Force (N)');
hold on
scatter(xn,yn);
disp(['Pour la fonction Fs, la valeur de R^2 est de: ' num2str(R) ' et la valeur RMS est de: ' num2str(RMS)]);

