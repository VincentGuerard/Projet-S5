%%
%Projet S5
%Aymen Mousli
% Linéraisation

%%
rabc = 95.2e-03;
bE1 = 13.029359254409743;
g = 9.81;
m = 0.442;
R = 3.6;
L = 0.115;

Xa = rabc;
Ya = 0;
Za = 0;
Xb = -rabc*sind(30);
Yb = rabc*cosd(30);
Zb = 0;
Xc = -rabc*sind(30);
Yc = -rabc*cosd(30);
Zc = 0;

aE0 = X1(1);
aE1 = X1(2);
aE2 = X1(3);
aE3 = X1(4);

as0 = X3(1);
as1 = X3(2);
as2 = X3(3);
as3 = X3(4);

Jp = 347e-06


%%

syms zo theta phi I Feq Deltaphi Deltatheta Deltazo deltaI

FEA = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xa*theta + Ya*phi) + aE2*(zo - Xa*theta + Ya*phi)^2 + aE3*(zo - Xa*theta + Ya*phi)^3);
FSA = -1/(as0 + as1*(zo - Xa*theta + Ya*phi) + as2*(zo - Xa*theta + Ya*phi)^2 + as3*(zo - Xa*theta + Ya*phi)^3);
FA = FEA+FSA;

alphaA(phi,theta) = diff(FA,phi,1);
betaA(phi,theta) = diff(FA,theta,1);
gammaA = diff(FA,zo,1);
sigmaA = diff(FA,I,1);

% DeltaFA = FAeq + alphaA(0,0) * Deltaphi + betaA(0,0)*Deltatheta + gammaA*Deltazo + sigmaA*deltaI

%%
%Condtions d'équilibre
syms iAeq iBeq iCeq zoeq xs ys

%Equilibre des forces

FCeq(xs,ys) = -m*g*((Xa*Yb + Xb*ys - xs*Yb - ys*Xa)/(Yb*Xa + Xb*Yc - Xc*Yb -Yc*Xa))
FAeq(xs,ys) = FCeq*((Xb*Yc)/(Xa*Yb) - Xc/Xa) + m*g*((Xb*ys)/(Xa*Yb) - (xs*Yb)/(Xa*Yb))
FBeq(xs,ys) = -( (Yc/Yb)*FCeq + m*g*(ys/Yb))

%Equilibre des courants

fA(zoeq) = -iAeq^2 - bE1*iAeq - ( FAeq(1,1)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))
fB(zoeq) = -iAeq^2 - bE1*iAeq - ( FBeq(1,1)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))
fC(zoeq) = -iAeq^2 - bE1*iAeq - ( FCeq(1,1)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))

Soluce_iAeq = vpasolve(fA(1)==0,iAeq)

%%
%Equilibre des tensions 

VAeq = R*iAeq
VBeq = R*iBeq
VCeq = R*iCeq

%Equilibre des angles

Phieq = 0;
Thetaeq =0;

%Equilibre des vitesse et accélération nulle

%%
% Equation linéaire

syms zo theta phi I Feq Deltaphi Deltatheta Deltazo deltaI

FEA = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xa*theta + Ya*phi) + aE2*(zo - Xa*theta + Ya*phi)^2 + aE3*(zo - Xa*theta + Ya*phi)^3);
FSA = -1/(as0 + as1*(zo - Xa*theta + Ya*phi) + as2*(zo - Xa*theta + Ya*phi)^2 + as3*(zo - Xa*theta + Ya*phi)^3);
FA = FEA+FSA;

alphaA(phi,theta) = diff(FA,phi,1);
betaA(phi,theta) = diff(FA,theta,1);
gammaA = diff(FA,zo,1);
sigmaA = diff(FA,I,1);

FEB = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xb*theta + Yb*phi) + aE2*(zo - Xb*theta + Yb*phi)^2 + aE3*(zo - Xb*theta + Yb*phi)^3);
FSB = -1/(as0 + as1*(zo - Xb*theta + Yb*phi) + as2*(zo - Xb*theta + Yb*phi)^2 + as3*(zo - Xb*theta + Yb*phi)^3);
FB = FEB+FSB;

alphaB(phi,theta) = diff(FB,phi,1);
betaB(phi,theta) = diff(FB,theta,1);
gammaB = diff(FB,zo,1);
sigmaB = diff(FB,I,1);

FEC = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xc*theta + Yc*phi) + aE2*(zo - Xc*theta + Yc*phi)^2 + aE3*(zo - Xc*theta + Yc*phi)^3);
FSC = -1/(as0 + as1*(zo - Xc*theta + Yc*phi) + as2*(zo - Xc*theta + Yc*phi)^2 + as3*(zo - Xc*theta + Yc*phi)^3);
FC = FEC+FSC;

alphaC(phi,theta) = diff(FC,phi,1);
betaC(phi,theta) = diff(FC,theta,1);
gammaC = diff(FC,zo,1);
sigmaC = diff(FC,I,1);

PP = [Yb/Jp*alphaB(0,0)+Yc/Jp*alphaC(0,0) Yb/Jp*betaB(0,0)+Yc/Jp*betaC(0,0) Yb/Jp*gamma(0,0)+Yc/Jp*gammaC(0,0) ]



























