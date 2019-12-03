%%
%Projet S5
%Aymen Mousli
% Linéraisation

%% Constantes
m = mP+0.040;
RR = 3.6;

aE0 = X1(1);
aE1 = X1(2);
aE2 = X1(3);
aE3 = X1(4);

as0 = X3(1);
as1 = X3(2);
as2 = X3(3);
as3 = X3(4);

%% ?????
syms zo theta phi I Feq Deltaphi Deltatheta Deltazo deltaI

FEA(I,theta,phi,zo) = (-I^2 + bE1*I)/(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3);
FSA(I,theta,phi,zo) = -1/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3);
FA(I,theta,phi,zo) = FEA+FSA;

alphaA(phi,theta) = diff(FA,phi,1);
betaA(phi,theta) = diff(FA,theta,1);
gammaA = diff(FA,zo,1);
sigmaA = diff(FA,I,1);

% DeltaFA = FAeq + alphaA(0,0) * Deltaphi + betaA(0,0)*Deltatheta + gammaA*Deltazo + sigmaA*deltaI

%% Condtions d'équilibre
syms iAeq iBeq iCeq zo xs ys

%Equilibre des forces
FCeq(xs,ys) = -m*g*((XA*YB + XB*ys - xs*YB - ys*XA)/(YB*XA + XB*YC - XC*YB -YC*XA));
FAeq(xs,ys) = FCeq(0.003,0.05)*((XB*YC)/(XA*YB) - XC/XA) + m*g*((XB*ys)/(XA*YB) - (xs*YB)/(XA*YB));
FBeq(xs,ys) = -( (YC/YB)*FCeq(0.003,0.05) + m*g*(ys/YB));

%Equilibre des courants
if FAeq(0.003,0.05) <0
    fA(theta,phi,zo) = -iAeq^2 + bE1*iAeq - ( FAeq(0.003,0.05)*(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3) + (aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3)/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3));
    
    Soluce_iAeq = vpasolve(fA(deg2rad(0.5),deg2rad(1),0.011)==0,iAeq);
    if Soluce_iAeq(1) <0 && imag(Soluce_iAeq(1))==0
        Soluce_iAeqfin = Soluce_iAeq(1);
    else
        Soluce_iAeqfin = Soluce_iAeq(2);
    end

end

if FAeq(0.003,0.05) >0
    fA(theta,phi,zo) = iAeq^2 + bE1*iAeq - ( FAeq(0.003,0.05)*(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3) + (aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3)/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3));
    
    Soluce_iAeq = vpasolve(fA(deg2rad(0.5),deg2rad(1),0.011)==0,iAeq);
    if Soluce_iAeq(1) >0 && imag(Soluce_iAeq(1))==0
        Soluce_iAeqfin = Soluce_iAeq(1);
    else
        Soluce_iAeqfin = Soluce_iAeq(2);
    end

end

if FBeq(0.003,0.05) <0
    fB(theta,phi,zo) = -iBeq^2 + bE1*iBeq - ( FBeq(0.003,0.05)*(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3) + (aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3)/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3));

    Soluce_iBeq = vpasolve(fB(deg2rad(0.5),deg2rad(1),0.011)==0,iBeq);
    if Soluce_iBeq(1) <0 && imag(Soluce_iBeq(1))==0
        Soluce_iBeqfin = Soluce_iBeq(1);
    else
        Soluce_iBeqfin = Soluce_iBeq(2);
    end

end

if FBeq(0.003,0.05) >0
    fB(theta,phi,zo) = iBeq^2 + bE1*iBeq - ( FBeq(0.003,0.05)*(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3) + (aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3)/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3));
    
    Soluce_iBeq = vpasolve(fB(deg2rad(0.5),deg2rad(1),0.011)==0,iBeq);
    if Soluce_iBeq(1) >0 && imag(Soluce_iBeq(1))==0
        Soluce_iBeqfin = Soluce_iBeq(1);
    else
        Soluce_iBeqfin = Soluce_iBeq(2);
    end

end

if FCeq(0.003,0.05) <0
    fC(theta,phi,zo) = -iCeq^2 + bE1*iCeq - ( FCeq(0.003,0.05)*(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3) + (aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3)/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3));

    Soluce_iCeq = vpasolve(fC(deg2rad(0.5),deg2rad(1),0.011)==0,iCeq);
    if Soluce_iCeq(1) <0 && imag(Soluce_iCeq(1))==0
        Soluce_iCeqfin = Soluce_iCeq(1);
    else
        Soluce_iCeqfin = Soluce_iCeq(2);
    end

end

if FCeq(0.003,0.05) >0
   fC(theta,phi,zo) = iCeq^2 + bE1*iCeq - ( FCeq(0.003,0.05)*(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3) + (aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3)/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3));

    Soluce_iCeq = vpasolve(fC(deg2rad(0.5),deg2rad(1),0.011)==0,iCeq);
    if Soluce_iCeq(1) >0 && imag(Soluce_iCeq(1))==0
        Soluce_iCeqfin = Soluce_iCeq(1);
    else
        Soluce_iCeqfin = Soluce_iCeq(2);
    end

end

%Équilibre des tensions 
VAeq = RR*Soluce_iAeqfin;
VBeq = RR*Soluce_iBeqfin;
VCeq = RR*Soluce_iCeqfin;

%Équilibre des angles
Phieq = 0;
Thetaeq =0;

%% Équation linéaire
syms zo theta phi Feq Deltaphi Deltatheta Deltazo deltaI Ia

FEA(phi,theta,I,zo) = (I^2 + bE1*I)/(aE0 + aE1*(zo - XA*theta + YA*phi) + aE2*(zo - XA*theta + YA*phi)^2 + aE3*(zo - XA*theta + YA*phi)^3);
FSA(phi,theta,I,zo) = -1/(as0 + as1*(zo - XA*theta + YA*phi) + as2*(zo - XA*theta + YA*phi)^2 + as3*(zo - XA*theta + YA*phi)^3);
FA(phi,theta,I,zo) = FEA(phi,theta,I,zo)+FSA(phi,theta,I,zo);

alphaA(phi,theta,I,zo) = diff(FA,phi,1);
betaA(phi,theta,I,zo) = diff(FA,theta,1);
gammaA(phi,theta,I,zo) = diff(FA,zo,1);
sigmaA(phi,theta,I,zo) = diff(FA,I,1);

FEB(phi,theta,I,zo) = (-I^2 + bE1*I)/(aE0 + aE1*(zo - XB*theta + YB*phi) + aE2*(zo - XB*theta + YB*phi)^2 + aE3*(zo - XB*theta + YB*phi)^3);
FSB(phi,theta,I,zo) = -1/(as0 + as1*(zo - XB*theta + YB*phi) + as2*(zo - XB*theta + YB*phi)^2 + as3*(zo - XB*theta + YB*phi)^3);
FB(phi,theta,I,zo) = FEB(phi,theta,I,zo)+FSB(phi,theta,I,zo);

alphaB(phi,theta,I,zo) = diff(FB,phi,1);
betaB(phi,theta,I,zo) = diff(FB,theta,1);
gammaB(phi,theta,I,zo) = diff(FB,zo,1);
sigmaB(phi,theta,I,zo) = diff(FB,I,1);

FEC(phi,theta,I,zo) = (I^2 + bE1*I)/(aE0 + aE1*(zo - XC*theta + YC*phi) + aE2*(zo - XC*theta + YC*phi)^2 + aE3*(zo - XC*theta + YC*phi)^3);
FSC(phi,theta,I,zo) = -1/(as0 + as1*(zo - XC*theta + YC*phi) + as2*(zo - XC*theta + YC*phi)^2 + as3*(zo - XC*theta + YC*phi)^3);
FC(phi,theta,I,zo) = FEC(phi,theta,I,zo)+FSC(phi,theta,I,zo);

alphaC(phi,theta,I,zo) = diff(FC,phi,1);
betaC(phi,theta,I,zo) = diff(FC,theta,1);
gammaC(phi,theta,I,zo) = diff(FC,zo,1);
sigmaC(phi,theta,I,zo) = diff(FC,I,1);

%% MATRICES
%Matrices de zéros et identitées
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

%Matrices de linéairisation
PP = [YB/Jx*alphaB(0,0,Soluce_iBeqfin,0.012)+YC/Jx*alphaC(0,0,Soluce_iCeqfin,0.012) YB/Jx*betaB(0,0,Soluce_iBeqfin,0.012)+YC/Jx*betaC(0,0,Soluce_iCeqfin,0.012) YB/Jx*gammaB(0,0,Soluce_iBeqfin,0.012)+YC/Jx*gammaC(0,0,Soluce_iCeqfin,0.012);
      -XA/Jy*alphaA(0,0,Soluce_iAeqfin,0.012)-XB/Jy*alphaB(0,0,Soluce_iBeqfin,0.012)-XC/Jy*alphaC(0,0,Soluce_iCeqfin,0.012) -XA/Jy*betaA(0,0,Soluce_iAeqfin,0.012)-XB/Jy*betaB(0,0,Soluce_iBeqfin,0.012)-XC/Jy*betaC(0,0,Soluce_iCeqfin,0.012) -XA/Jy*gammaA(0,0,Soluce_iAeqfin,0.012)-XB/Jy*gammaB(0,0,Soluce_iBeqfin,0.012)-XC/Jy*gammaC(0,0,Soluce_iCeqfin,0.012);
      alphaA(0,0,Soluce_iAeqfin,0.012)/m+alphaB(0,0,Soluce_iBeqfin,0.012)/m+alphaC(0,0,Soluce_iCeqfin,0.012)/m betaA(0,0,Soluce_iAeqfin,0.012)/m+betaB(0,0,Soluce_iBeqfin,0.012)/m+betaC(0,0,Soluce_iCeqfin,0.012)/m gammaA(0,0,Soluce_iAeqfin,0.012)/m+gammaB(0,0,Soluce_iBeqfin,0.012)/m+gammaC(0,0,Soluce_iCeqfin,0.012)/m];
PS = [0 mS*g/Jx;
      -mS*g/Jy 0;
      0 0];
PC = [0 YB/Jx*sigmaB(0,0,Soluce_iBeqfin,0.012) YC/Jx*sigmaC(0,0,Soluce_iCeqfin,0.012);
      -XA/Jx*sigmaA(0,0,Soluce_iAeqfin,0.012) -XB/Jx*sigmaB(0,0,Soluce_iBeqfin,0.012) -XC/Jx*sigmaC(0,0,Soluce_iCeqfin,0.012);
      sigmaA(0,0,Soluce_iAeqfin,0.012)/m sigmaB(0,0,Soluce_iBeqfin,0.012)/m sigmaC(0,0,Soluce_iCeqfin,0.012)/m];
SP = [0 -5/7*g 0;
      5/7*g 0 0];
CC = [-RR/LL 0 0;
      0 -RR/LL 0;
      0 0 -RR/LL];
CV = [1/LL 0 0;
      0 1/LL 0;
      0 0 1/LL];
Tdef = [YD -XD 1;
        YE -XE 1;
        YF -XF 1];

%Matrices A,B,C et D du modèle linéaire
A = [mat_Zero33 mat_One33 mat_Zero32 mat_Zero32 mat_Zero33;
     PP mat_Zero33 PS mat_Zero32 PC;
     mat_Zero23 mat_Zero23 mat_Zero22 mat_One22 mat_Zero23;
     SP mat_Zero23 mat_Zero22 mat_Zero22 mat_Zero23;
     mat_Zero33 mat_Zero33 mat_Zero32 mat_Zero32 CC];

B = [mat_Zero33; mat_Zero33; mat_Zero23; mat_Zero23; CV];

C = [Tdef mat_Zero33 mat_Zero34 mat_Zero33;
     mat_Zero43 mat_Zero43 mat_One44 mat_Zero43];
 
D = mat_Zero73;

A = double(A);
B = double(B);
C = double(C);
D = double(D);

%%
ZZABC = TABC'*[0;0;0.012]
ZZDEF = TDEF'*[deg2rad(1);deg2rad(0.5);0.011]