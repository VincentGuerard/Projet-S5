% Trajectoire
function [Pi, Ltr, E, Vr, Traj, tt] = trajectoire(position, vAB, Ts)

xn = position(:,1);
yn = position(:,2);

%% -------------------------- Interpolation -------------------------------

% Coefficients du Polyn�me d'interpolation (f)
for i = 1:length(xn)    
    P(:,i) =  xn.^(i-1);
end

A = pinv(P)*yn;

for j = 1:length(xn)
    Pi(1,j) = A(j);
end

Pi = fliplr(Pi);

% Coefficients du Polyn�me d'interpolation (f')
N = length(Pi);
for i=1 : N-1
    Pi_der(i) = Pi(i)*(N-i);
end

% Coefficients du Polyn�me d'interpolation (f'')
Nr = length(Pi_der);
for i=1 : Nr-1
    Pi_der2(i) = Pi_der(i)*(Nr-i);
end

%% -------------------- Longueur de la Trajectoire ------------------------
M = 101;
xrow = sortrows(xn);
xm = xrow(1): xrow(end)/M : xrow(end);
h = xrow(end)/M;

y = polyval(Pi,xm);           % ----> f
yp = polyval(Pi_der,xm);      % ----> f'
ypp = polyval(Pi_der2,xm);    % ----> f''
g = sqrt(1+yp.^2);            % ----> g
gp = ypp.*yp./g;              % ----> g'

for i = 2 : length(g)
    Longueur(i) =(h/2)*(g(1) + g(i) + 2*sum(g(2:i-1)));
end

% Variable de Sortie
Ltr(:,1) = xm;
Ltr(:,2) = Longueur;

%% ------------------------- Erreur d'int�gration --------------------------

gpa = gp(1);                  % ----> g'(a)
gpb = gp(end);                % ----> g'(b)

E = (h^2/12)*(gpb-gpa);

%% ---------------------- �chantilonnage des points -----------------------

% Calcul du nombre de points O et de la vitesse r�elle
dt = vAB * Ts;
O = round(Ltr(end,2)/dt);
% O = 10;
dt = Ltr(end,2)/O;
Vr = dt/Ts;

% Pr�cision du Newton-Raphson
tol = 1e-08;
iteMax = 101;

% Conditions initiales
Traj = zeros(O,2);
D = Longueur(end)/(O-1);
a0 = 0;
a = a0;

for i= 1 : O
    b = a + 0.001;
    y = polyval(Pi,a);
    yp = polyval(Pi_der,a);
    g = sqrt(1+yp.^2);
    h = 0.1;
    fn = 0 - D;
    fpn = g;
    ite = 0;
    while abs(fn) > tol && ite < iteMax
        b = b - fn/fpn;
        h = (b-a)/(M-1);
        x = a:h:b;
        y = polyval(Pi,x);
        yp = polyval(Pi_der,x);
        g = sqrt(1+yp.^2);
        fn = (h/2)*(g(1) + g(end) + 2*sum(g(2:end-1)))- D;
        fpn = g(end);
        ite = ite+1;
    end    
    Traj(i,1) = b;
    Traj(i,2) = polyval(Pi,b);
    
    a = b;
end

% Calcul du temps total
tt = Ltr(end,2)/Vr;

end