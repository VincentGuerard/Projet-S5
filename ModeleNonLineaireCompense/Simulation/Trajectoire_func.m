% Trajectoire


function [Pi, Ltr, E, Vr, Traj, tt] = Trajectoire_func(position, vAB, Ts)

xn = position(:,1);
yn = position(:,2);

%% -------------------------- Interpolation -------------------------------

% Coefficients du Polynôme d'interpolation (f)
for i = 1:length(xn)    
    P(:,i) =  xn.^(i-1);
end
A = pinv(P)*yn;

for j = 1:length(xn)
    Pi(1,j) = A(j);
end

Pi = fliplr(Pi);

% Coefficients du Polynôme d'interpolation (f')
N = length(Pi);
for i=1 : N-1
    Pi_der(i) = Pi(i)*(N-i);
end

% Coefficients du Polynôme d'interpolation (f'')
Nr = length(Pi_der);
for i=1 : Nr-1
    Pi_der2(i) = Pi_der(i)*(Nr-i);
end

%% -------------------- Longueur de la Trajectoire ------------------------
M = 101;
xrow = sortrows(xn);
h = (xrow(end)-xrow(1))/(M-1);
xm = xrow(1): h : xrow(end);

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

%% ------------------------- Erreur d'intégration --------------------------

gpa = gp(1);                  % ----> g'(a)
gpb = gp(end);                % ----> g'(b)

E = (h^2/12)*(gpb-gpa);

%% ---------------------- Échantilonnage des points -----------------------

% Calcul du nombre de points O
dt = vAB * Ts;
O = round(Ltr(end,2)/dt);

% Calcul de la vitesse réelle
dt = Ltr(end,2)/O;
Vr = dt/Ts;

% Calcul du temps total
tt = Ltr(end,2)/Vr;

% Précision du Newton-Raphson
tol = 1e-08;
iteMax = 101;

% Conditions initiales
Traj = zeros(O,2);
D = Longueur(end)/(O-1);
a0 = min(xn);
a = a0;

% Calcul de la prochaine borne d'intégration

    b = a + 0.001;               % ----> borne b légèrement supérieure
for i= 1 : O  
    yp = polyval(Pi_der,a);       % ----> f'
    g = sqrt(1+yp.^2);            % ----> g
    fn = -D;                      % ----> g(0) - D
    fpn = g;                      % ----> dL/db = g(b)
    ite = 0;
    while abs(fn) > tol && ite < iteMax
        b = b - fn/fpn;
        h = (b-a)/(M-1);
        x = a:h:b;
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

end