%% Mod�le Decoupl�
clear all
close all
clc

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

syms Fek Fsk alpha_a alpha_b alpha_c beta_a beta_b beta_c gamma_a gamma_b gamma_c sigma_a sigma_b sigma_c YA YB YC XA XB XC m Jpx Jpy masseS g
syms Jsx Jsy Ra Rb Rc La Lb Lc YD YE YF XD XE XF

PP = [YB/Jpx*alpha_b+YC/Jpx*alpha_c YB/Jpx*beta_b+YC/Jpx*beta_c YB/Jpx*gamma_b+YC/Jpx*gamma_c;
      -XA/Jpy*alpha_a-XB/Jpy*alpha_b-XC/Jpy*alpha_c XA/Jpy*beta_a-XB/Jpy*beta_b-XC/Jpy*beta_c -XA/Jpy*gamma_a-XB/Jpy*gamma_b-XC/Jpy*gamma_c;
      alpha_a/m+alpha_b/m+alpha_c/m beta_a/m+beta_b/m+beta_c/m gamma_a/m+gamma_b/m+gamma_c/m];
PS = [0 masseS*g/Jsx;
      -masseS*g/Jsy 0;
      0 0];
PC = [0 YB/Jpx*sigma_b YC/Jpx*sigma_c;
      -XA/Jpx*sigma_a -XB/Jpx*sigma_b -XC/Jpx*sigma_c;
      sigma_a/m sigma_b/m sigma_c/m];
SP = [0 -5/7*g 0;
      5/7*g 0 0];
CC = [-Ra/La 0 0;
      0 -Rb/Lb 0;
      0 0 -Rc/Lc];
CV = [1/La 0 0;
      0 1/Lb 0;
      0 0 1/Lc];
Tdef = [YD -XD 1;
        YE -XE 1;
        YF -XF 1];
Tabc = [YA YB YC;
        -XA -XB -XC;
        1 1 1];

A = [mat_Zero33 mat_One33 mat_Zero32 mat_Zero32 mat_Zero33;
     PP mat_Zero33 PS mat_Zero32 PC;
     mat_Zero23 mat_Zero23 mat_Zero22 mat_One22 mat_Zero23;
     SP mat_Zero23 mat_Zero22 mat_Zero22 mat_Zero23;
     mat_Zero33 mat_Zero33 mat_Zero32 mat_Zero32 CC];

B = [mat_Zero33; mat_Zero33; mat_Zero23; mat_Zero23; CV];
C = [Tdef mat_Zero33 mat_Zero34 mat_Zero33;
     mat_Zero43 mat_Zero43 mat_One44 mat_Zero43];
D = [mat_Zero73];
    
%% Remplacement des variables
Rs = 3.9/1000;
masseS = 8/1000;
g = 9.81;
Js = 2*masseS*Rs^2/5;

Ra = 1;
Rb = 1;
Rc = 1;
La = 1;
Lb = 1;
Lc = 1;

r_abc = 95.2/1000;
XA = r_abc;
XB = -r_abc*sind(30);
XC = -r_abc*sind(30);
XD = 1;
XE = 1;
XF = 1;
YA = 0;
YB = r_abc*cosd(30);
YC = -r_abc*cosd(30);
YD = 1;
YE = 1;
YF = 1;
Jxy = 1347*10^-6;
Jpx = Jxy;
Jpy = Jxy;
Jsx = 1;
Jsy = 1;
masseP = 442/1000
m = masseP;

be1 = 1;
ae0 = 1;
ae1 = 1;
ae2 = 1;
ae3 = 1;
as0 = 1;
as1 = 1;
as2 = 1;
as3 = 1;

L = 1;
R = 1;
V = 1;

alpha_a = 1;
alpha_b = 1;
alpha_c = 1;
beta_a = 1;
beta_b = 1;
beta_c = 1;
gamma_a = 1;
gamma_b = 1;
gamma_c = 1;
sigma_a = 1;
sigma_b = 1;
sigma_c = 1;

A = double(subs(A))
B = double(subs(B))
C = double(subs(C))

%% Simulation 
u = [1; 1; 1];

%% Affichage des r�sultats



