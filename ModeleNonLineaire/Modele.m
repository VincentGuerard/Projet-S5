%% Modelisation de la partie électrique des actionneurs
syms ia ib ic Ra Rb Rc La Lb Lc Va Vb Vc s

A = [-Ra/La 0 0;
     0 -Rb/Lb 0;
     0 0 -Rc/Lc];
B = [1/La 0 0;
     0 1/Lb 0;
     0 0 1/Lc];
C = [1 1 1];
D = 0
I = [s 0 0;
     0 s 0;
     0 0 s];
FT = C*(inv(I-A)*B)+D

[numA, denA] = numden(FT(1))

L = 1
R =1
V = 1
%% Modelisation de la partie electromécanique des actionneurs
syms Fek Fsk Fk ik be1 ae0 ae1 ae2 ae3 as0 as1 as2 as3 zk
be1 = 1;
ae0 = 1;
ae1 = 1;
ae2 = 1;
ae3 = 1;
as0 = 1;
as1 = 1;
as2 = 1;
as3 = 1;


%% Dynamique de la plaque
r_abc = 95.2/1000;
XA = r_abc;
YA = 0;
XB = -r_abc*sind(30);
YB = r_abc*cosd(30);
XC = -r_abc*sind(30);
YC = -r_abc*cosd(30);
Jxy = 1347*10^-6;
masseP = 442/1000


%% Dynamique de la sphère
Rs = 3.9/1000;
masseS = 8/1000;
g = 9.81;
Js = 2*masseS*Rs^2/5;
