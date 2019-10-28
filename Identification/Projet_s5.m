%%
% Projet de S5
% Recherche des coefficient pour les forces fk
% Aymen Mousli
% 23/10/2019 11:45 v1

%%
% clear all
% close all
% clc

%%
%Création des variables
be1 = 13.029359254409743;

% Fs = load('Fs')
% Fe = load('Fe_attraction')
Fe = load("..\DonneesExperimentales\Fe_attraction.mat")
Fs = load("..\DonneesExperimentales\Fs.mat")
bE1 = 13.029359254409743;
Offset = 0:0.1:15;
I1 = -(ones(size(Fe.Fe_m1A)));
I2 = -(2*ones(size(Fe.Fe_m2A)));
Y1 = -(I1.^2 - bE1*I1);
Y2 = -(I2.^2 - bE1*I2);
%%
% Pour I = -1
for i=1 : length(Offset)
    N = length(Fe.Fe_m1A);
    FeOffset1 = Fe.Fe_m1A + Offset(i);
    zOffset1 = Fe.z_m1A + Offset(i);
    A = [FeOffset1 FeOffset1.*zOffset1 (FeOffset1.*(zOffset1.^2)) (FeOffset1.*(zOffset1.^3))];
    X = pinv(A)*Y1;
    ybarre = mean(FeOffset1);
    yappro = Y1./(X(1)+(X(2).*zOffset1)+(X(3).*zOffset1.*zOffset1)+(X(4).*zOffset1.*zOffset1.*zOffset1));
    RR(i)= 1-(sum((FeOffset1-yappro).^2)/sum((FeOffset1-ybarre).^2));
end

%%
% Recherche meilleur offset
max = RR(1)
for i=1 : length(RR)
    if RR(i)>max
        max = RR(i);
        position = i;
    end
    
end
%%
%Appro

A = [Fe.Fe_m1A Fe.Fe_m1A.*Fe.z_m1A (Fe.Fe_m1A.*(Fe.z_m1A.^2)) (Fe.Fe_m1A.*(Fe.z_m1A.^3))]
X1 = pinv(A)*Y1
yappro = Y1./(X1(1)+(X1(2).*Fe.z_m1A)+(X1(3).*Fe.z_m1A.*Fe.z_m1A)+(X1(4).*Fe.z_m1A.*Fe.z_m1A.*Fe.z_m1A))

%%
figure
plot(Fe.z_m1A,Fe.Fe_m1A)
hold on 
plot(Fe.z_m1A,yappro)

%%
% Pour I=-2

for i=1 : length(Offset)
    N = length(Fe.Fe_m2A);
    FeOffset2 = Fe.Fe_m2A + Offset(i);
    zOffset2 = Fe.z_m2A + Offset(i);
    A = [FeOffset2 FeOffset2.*zOffset2 (FeOffset2.*(zOffset2.^2)) (FeOffset2.*(zOffset2.^3))];
    X = pinv(A)*Y2;
    ybarre = mean(FeOffset2);
    yappro = Y2./(X(1)+(X(2).*zOffset2)+(X(3).*zOffset2.*zOffset2)+(X(4).*zOffset2.*zOffset2.*zOffset2));
    RR(i)= 1-(sum((FeOffset2-yappro).^2)/sum((FeOffset2-ybarre).^2));
end

max = RR(1)
for i=1 : length(RR)
    if RR(i)>max
        max = RR(i);
        position = i;
    end
    
end

%%
%Appro


A = [Fe.Fe_m2A Fe.Fe_m2A.*Fe.z_m2A (Fe.Fe_m2A.*(Fe.z_m2A.^2)) (Fe.Fe_m2A.*(Fe.z_m2A.^3))];
X2 = pinv(A)*Y2;
yappro = Y2./(X2(1)+(X2(2).*Fe.z_m2A)+(X2(3).*Fe.z_m2A.*Fe.z_m2A)+(X2(4).*Fe.z_m2A.*Fe.z_m2A.*Fe.z_m2A));

figure
plot(Fe.z_m2A,Fe.Fe_m2A)
hold on 
plot(Fe.z_m2A,yappro)

%%
% Pour FS
Y3 = -(ones(size(Fs.Fs)));
A = [Fs.Fs Fs.Fs.*Fs.z_pos (Fs.Fs.*(Fs.z_pos.^2)) (Fs.Fs.*(Fs.z_pos.^3))];
X3 = pinv(A)*Y3;
yappro = Y3./(X3(1)+(X3(2).*Fs.z_pos)+(X3(3).*Fs.z_pos.*Fs.z_pos)+(X3(4).*Fs.z_pos.*Fs.z_pos.*Fs.z_pos));

figure
plot(Fs.z_pos,Fs.Fs)
hold on 
plot(Fs.z_pos,yappro)

%% Mapping vers variables du simulink
as0 = X3(1);
as1 = X3(2);
as2 = X3(3);
as3 = X3(4);

ae0 = X1(1);
ae1 = X1(2);
ae2 = X1(3);
ae3 = X1(4);

