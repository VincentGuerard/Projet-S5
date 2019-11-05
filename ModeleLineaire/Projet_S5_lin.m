%%
%Projet S5
%Aymen Mousli
% Linéraisation

%%
Rs = 3.9/1000;
masseS = 8/1000;
g = 9.81;
Js = 2*masseS*Rs^2/5;

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

XA = rabc;
YA = 0;
ZA = 0;
XB = -rabc*sind(30);
YB = rabc*cosd(30);
ZB = 0;
XC = -rabc*sind(30);
YC = -rabc*cosd(30);
ZC = 0;

rdef = 80/1000;
XD = rdef*sind(30);
XE = -rdef;
XF = rdef*sind(30);
YD = rdef*cosd(30);
YE = 0;
YF = -rdef*cosd(30);

aE0 = X1(1);
aE1 = X1(2);
aE2 = X1(3);
aE3 = X1(4);

as0 = X3(1);
as1 = X3(2);
as2 = X3(3);
as3 = X3(4);

Jpx = 347e-06

Jpy = 347e-06


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
FAeq(xs,ys) = FCeq(0,0)*((Xb*Yc)/(Xa*Yb) - Xc/Xa) + m*g*((Xb*ys)/(Xa*Yb) - (xs*Yb)/(Xa*Yb))
FBeq(xs,ys) = -( (Yc/Yb)*FCeq(0,0) + m*g*(ys/Yb))

%Equilibre des courants

if FAeq(0,0) <0
    fA(zoeq) = -iAeq^2 + bE1*iAeq - ( FAeq(0,0)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))
    
    Soluce_iAeq = vpasolve(fA(0.015)==0,iAeq)
    if Soluce_iAeq(1) <0 && imag(Soluce_iAeq(1))==0
        Soluce_iAeqfin = Soluce_iAeq(1)
    else
        Soluce_iAeqfin = Soluce_iAeq(2)
    end

end

if FAeq(0,0) >0
    fA(zoeq) = iAeq^2 + bE1*iAeq - ( FAeq(0,0)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))
    
    Soluce_iAeq = vpasolve(fA(0.015)==0,iAeq)
    if Soluce_iAeq(1) >0 && imag(Soluce_iAeq(1))==0
        Soluce_iAeqfin = Soluce_iAeq(1)
    else
        Soluce_iAeqfin = Soluce_iAeq(2)
    end

end

if FBeq(0,0) <0
    fB(zoeq) = -iBeq^2 + bE1*iBeq - ( FBeq(0,0)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))

    Soluce_iBeq = vpasolve(fB(0.015)==0,iBeq)
    if Soluce_iBeq(1) <0 && imag(Soluce_iBeq(1))==0
        Soluce_iBeqfin = Soluce_iBeq(1)
    else
        Soluce_iBeqfin = Soluce_iBeq(2)
    end

end

if FBeq(0,0) >0
    fB(zoeq) = iBeq^2 + bE1*iBeq - ( FBeq(0,0)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))
    
    Soluce_iBeq = vpasolve(fB(0.015)==0,iBeq)
    if Soluce_iBeq(1) >0 && imag(Soluce_iBeq(1))==0
        Soluce_iBeqfin = Soluce_iBeq(1)
    else
        Soluce_iBeqfin = Soluce_iBeq(2)
    end

end

if FCeq(0,0) <0
    fC(zoeq) = -iCeq^2 + bE1*iCeq - ( FCeq(0,0)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))

    Soluce_iCeq = vpasolve(fC(0.015)==0,iCeq)
    if Soluce_iCeq(1) <0 && imag(Soluce_iCeq(1))==0
        Soluce_iCeqfin = Soluce_iCeq(1)
    else
        Soluce_iCeqfin = Soluce_iCeq(2)
    end

end

if FCeq(0,0) >0
   fC(zoeq) = iCeq^2 + bE1*iCeq - ( FCeq(0,0)*(aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3) + (aE0 + aE1*zoeq + aE2*zoeq^2 + aE3*zoeq^3)/(as0 + as1*zoeq + as2*zoeq^2 + as3*zoeq^3))

    Soluce_iCeq = vpasolve(fC(0.015)==0,iCeq)
    if Soluce_iCeq(1) >0 && imag(Soluce_iCeq(1))==0
        Soluce_iCeqfin = Soluce_iCeq(1)
    else
        Soluce_iCeqfin = Soluce_iCeq(2)
    end

end

%%
%Equilibre des tensions 

VAeq = R*Soluce_iAeqfin
VBeq = R*Soluce_iBeqfin
VCeq = R*Soluce_iCeqfin

%Equilibre des angles

Phieq = 0;
Thetaeq =0;

%Equilibre des vitesse et accélération nulle

%%
% Equation linéaire

syms zo theta phi Feq Deltaphi Deltatheta Deltazo deltaI Ia

FEA(phi,theta,I,zo) = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xa*theta + Ya*phi) + aE2*(zo - Xa*theta + Ya*phi)^2 + aE3*(zo - Xa*theta + Ya*phi)^3);
FSA(phi,theta,I,zo) = -1/(as0 + as1*(zo - Xa*theta + Ya*phi) + as2*(zo - Xa*theta + Ya*phi)^2 + as3*(zo - Xa*theta + Ya*phi)^3);
FA(phi,theta,I,zo) = FEA(phi,theta,I,zo)+FSA(phi,theta,I,zo);

alphaA(phi,theta,I,zo) = diff(FA,phi,1);
betaA(phi,theta,I,zo) = diff(FA,theta,1);
gammaA(phi,theta,I,zo) = diff(FA,zo,1);
sigmaA(phi,theta,I,zo) = diff(FA,I,1);

FEB(phi,theta,I,zo) = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xb*theta + Yb*phi) + aE2*(zo - Xb*theta + Yb*phi)^2 + aE3*(zo - Xb*theta + Yb*phi)^3);
FSB(phi,theta,I,zo) = -1/(as0 + as1*(zo - Xb*theta + Yb*phi) + as2*(zo - Xb*theta + Yb*phi)^2 + as3*(zo - Xb*theta + Yb*phi)^3);
FB(phi,theta,I,zo) = FEB(phi,theta,I,zo)+FSB(phi,theta,I,zo);

alphaB(phi,theta,I,zo) = diff(FB,phi,1);
betaB(phi,theta,I,zo) = diff(FB,theta,1);
gammaB(phi,theta,I,zo) = diff(FB,zo,1);
sigmaB(phi,theta,I,zo) = diff(FB,I,1);

FEC(phi,theta,I,zo) = (-I^2 + bE1*I)/(aE0 + aE1*(zo - Xc*theta + Yc*phi) + aE2*(zo - Xc*theta + Yc*phi)^2 + aE3*(zo - Xc*theta + Yc*phi)^3);
FSC(phi,theta,I,zo) = -1/(as0 + as1*(zo - Xc*theta + Yc*phi) + as2*(zo - Xc*theta + Yc*phi)^2 + as3*(zo - Xc*theta + Yc*phi)^3);
FC(phi,theta,I,zo) = FEC(phi,theta,I,zo)+FSC(phi,theta,I,zo);

alphaC(phi,theta,I,zo) = diff(FC,phi,1);
betaC(phi,theta,I,zo) = diff(FC,theta,1);
gammaC(phi,theta,I,zo) = diff(FC,zo,1);
sigmaC(phi,theta,I,zo) = diff(FC,I,1);

%%

mat_One33 = eye(3);
mat_One44 = eye(4);
mat_One22 = eye(2);
mat_Zero33 = zeros(3);
mat_Zero32 = zeros(3,2);
mat_Zero23 = zeros(2,3);
mat_Zero22 = zeros(2,2);
mat_Zero34 = zeros(3,4);
mat_Zero43 = zeros(4,3);
mat_Zero73 = zeros(7,3);


PP = [YB/Jpx*alphaB(0,0,Soluce_iBeqfin,0.015)+YC/Jpx*alphaC(0,0,Soluce_iCeqfin,0.015) YB/Jpx*betaB(0,0,Soluce_iBeqfin,0.015)+YC/Jpx*betaC(0,0,Soluce_iCeqfin,0.015) YB/Jpx*gammaB(0,0,Soluce_iBeqfin,0.015)+YC/Jpx*gammaC(0,0,Soluce_iCeqfin,0.015);
      -XA/Jpy*alphaA(0,0,Soluce_iAeqfin,0.015)-XB/Jpy*alphaB(0,0,Soluce_iBeqfin,0.015)-XC/Jpy*alphaC(0,0,Soluce_iCeqfin,0.015) -XA/Jpy*betaA(0,0,Soluce_iAeqfin,0.015)-XB/Jpy*betaB(0,0,Soluce_iBeqfin,0.015)-XC/Jpy*betaC(0,0,Soluce_iCeqfin,0.015) -XA/Jpy*gammaA(0,0,Soluce_iAeqfin,0.015)-XB/Jpy*gammaB(0,0,Soluce_iBeqfin,0.015)-XC/Jpy*gammaC(0,0,Soluce_iCeqfin,0.015);
      alphaA(0,0,Soluce_iAeqfin,0.015)/m+alphaB(0,0,Soluce_iBeqfin,0.015)/m+alphaC(0,0,Soluce_iCeqfin,0.015)/m betaA(0,0,Soluce_iAeqfin,0.015)/m+betaB(0,0,Soluce_iBeqfin,0.015)/m+betaC(0,0,Soluce_iCeqfin,0.015)/m gammaA(0,0,Soluce_iAeqfin,0.015)/m+gammaB(0,0,Soluce_iBeqfin,0.015)/m+gammaC(0,0,Soluce_iCeqfin,0.015)/m];
PS = [0 masseS*g/Jpx;
      -masseS*g/Jpy 0;
      0 0];
PC = [0 YB/Jpx*sigmaB(0,0,Soluce_iBeqfin,0.015) YC/Jpx*sigmaC(0,0,Soluce_iCeqfin,0.015);
      -XA/Jpx*sigmaA(0,0,Soluce_iAeqfin,0.015) -XB/Jpx*sigmaB(0,0,Soluce_iBeqfin,0.015) -XC/Jpx*sigmaC(0,0,Soluce_iCeqfin,0.015);
      sigmaA(0,0,Soluce_iAeqfin,0.015)/m sigmaB(0,0,Soluce_iBeqfin,0.015)/m sigmaC(0,0,Soluce_iCeqfin,0.015)/m];
SP = [0 -5/7*g 0;
      5/7*g 0 0];
CC = [-R/L 0 0;
      0 -R/L 0;
      0 0 -R/L];
CV = [1/L 0 0;
      0 1/L 0;
      0 0 1/L];
Tdef = [YD -XD 1;
        YE -XE 1;
        YF -XF 1];

A = [mat_Zero33 mat_One33 mat_Zero32 mat_Zero32 mat_Zero33;
     PP mat_Zero33 PS mat_Zero32 PC;
     mat_Zero23 mat_Zero23 mat_Zero22 mat_One22 mat_Zero23;
     SP mat_Zero23 mat_Zero22 mat_Zero22 mat_Zero23;
     mat_Zero33 mat_Zero33 mat_Zero32 mat_Zero32 CC];

B = [mat_Zero33; mat_Zero33; mat_Zero23; mat_Zero23; CV];
C = [Tdef mat_Zero33 mat_Zero34 mat_Zero33;
     mat_Zero43 mat_Zero43 mat_One44 mat_Zero43];
D = [mat_Zero73];

A = double(A);
B = double(B);
C = double(C);
D = double(D);



















