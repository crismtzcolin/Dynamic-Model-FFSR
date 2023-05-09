% Robot manipulador espacial de 5gdl:
syms PBx PBy q0 q1 q2 
syms DPBx DPBy Dq0 Dq1 Dq2 
syms DDPBx DDPBy DDq0 DDq1 DDq2 
syms m0 m1 m2 Pb L1 L2 a1 a2 I0xx I0yy I0zz I1xx I1yy I1zz I2xx I2yy I2zz tau0 tau1 tau2 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Matrices de inercia   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I0=[I0xx, 0, 0;
    0, I0yy, 0;
    0, 0, I0zz];

I1=[I1xx, 0, 0;
    0, I1yy, 0;
    0, 0, I1zz];

I2=[I2xx, 0, 0;
    0, I2yy, 0;
    0, 0, I2zz];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Vectores de velocidad lineal  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v0=[DPBx;
    DPBy;
    0];
v1=[DPBx-a1*sin(q0+q1)*(Dq0+Dq1);
    DPBy+a1*cos(q0+q1)*(Dq0+Dq1);
    0];
v2=[DPBx-L1*sin(q0+q1)*(Dq0+Dq1)-a2*sin(q0+q1+q2)*(Dq0+Dq1+Dq2);
    DPBy+L1*cos(q0+q1)*(Dq0+Dq1)+a2*cos(q0+q1+q2)*(Dq0+Dq1+Dq2);
    0];

% v0=[DPBx;
%     DPBy;
%     0];
% v1=[-a1*sin(q0+q1)*(Dq0+Dq1);
%     a1*cos(q0+q1)*(Dq0+Dq1);
%     0];
% v2=[-a2*sin(q0+q1+q2)*(Dq0+Dq1+Dq2);
%     a2*cos(q0+q1+q2)*(Dq0+Dq1+Dq2);
%     0];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Vectores de velocidad angular  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0=[0;
    0;
    Dq0];
w1=[0;
    0;
    Dq1];
w2=[0;
    0;
    Dq2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Declaración de variables adicionales  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

O3x3 = [0 , 0 , 0;
    0 , 0 , 0;
    0 , 0 , 0];

E = [1 , 0 , 0;
    0 , 1 , 0;
    0 , 0 , 1];

wc=m0+m1+m2;

dotPhi= [Dq1;
    Dq2];

SRb = [0, 0, PBy;
    0, 0, -PBx;
    -PBy,PBx,0];

SRc = [0, 0, (PBy+L1*sin(q0+q1)+a2*sin(q0+q1+q2));
    0, 0, -(PBx+L1*cos(q0+q1)+a2*cos(q0+q1+q2));
    -(PBx+L1*sin(q0+q1)+a2*sin(q0+q1+q2)),(PBy+L1*cos(q0+q1)+a2*cos(q0+q1+q2)),0];

Srr=[0, 0, (PBy+L1*sin(q0+q1)+a2*sin(q0+q1+q2));
    0, 0, -(PBx+L1*cos(q0+q1)+a2*cos(q0+q1+q2));
    -(PBy+L1*sin(q0+q1)+a2*sin(q0+q1+q2)),(PBx+L1*cos(q0+q1)+a2*cos(q0+q1+q2)),0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Jacobiano lineal del RM wrt B  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jv1=[-a2*sin(q0+q1+q2)-L1*sin(q0+q1);
    a2*cos(q0+q1+q2)+L1*cos(q0+q1);
    0];
Jv2=[-a2*sin(q0+q1+q2);
    a2*cos(q0+q1+q2);
    0];

k0=[0;
    0;
    1];

k1=[0;
    0;
    1];

k2=[0;
    0;
    1];


mJv=[m1*Jv1, m2*Jv2];
mJw=[I1*k1+m1*SRc*Jv1,I2*k2+m2*SRc*Jv2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Matriz de transformación N  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = -inv([wc*E, -wc*Srr; wc*SRc, I0+wc*SRc*Srr])*[[mJv;mJw]];
disp('Matriz N')
N=simplify(N);
disp(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   LAGRANGIANO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T   = ([N*dotPhi].')*[m0*E , O3x3;O3x3 , I0]*[N*dotPhi]+(v1.'*m1*v1+w1.'*m1*w1)+(v2.'*m2*v2+w2.'*m2*w2);
T=simplify(T)
disp('Energía cinética')
disp(T)

disp('Lagrangiano')

V   = 0;
L   = T - V;
disp(L)

fid = fopen('L.txt', 'wt');
fprintf(fid, '%s\n', char(L));
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Modelado Dinámico por Euler Lagrange  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Fase 1')
q  = [q1, q2];
Dq = [Dq1, Dq2];
tt = linspace(0,5,500);
disp('Fase 2')
Eq = LagrangeDynamicEqDeriver(L, q, Dq);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Fase 3')
disp("Modelo dinámico:")
disp(Eq)

fid = fopen('Modelo.txt', 'wt');
fprintf(fid, '%s\n', char(Eq));
fclose(fid);

fid = fopen('Eq1.txt', 'wt');
fprintf(fid, '%s\n', char(Eq(1)));
fclose(fid);

fid = fopen('Eq2.txt', 'wt');
fprintf(fid, '%s\n', char(Eq(2)));
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Fase 4')
disp("Matriz de inercia: ")
M=[diff(Eq(1),DDq1), diff(Eq(1),DDq2);
   diff(Eq(2),DDq1), diff(Eq(2),DDq2)];
disp(M)

fid = fopen('M.txt', 'wt');
fprintf(fid, '%s\n', char(M));
fclose(fid);

fid = fopen('M Python.txt', 'wt');
fprintf(fid, '%s\n', char(M));
fclose(fid);

fid = fopen('Inercia.m', 'wt');
fprintf(fid, '%s\n', char(M));
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Fase 5')
disp("Matriz de coriolis: ")
C=[diff(Eq(1),Dq1),diff(Eq(1),Dq2);
   diff(Eq(2),Dq1),diff(Eq(2),Dq2)];
disp(C)

fid = fopen('C.txt', 'wt');
fprintf(fid, '%s\n', char(C));
fclose(fid);

fid = fopen('C Python.txt', 'wt');
fprintf(fid, '%s\n', char(C));
fclose(fid);

fid = fopen('Coriolis.m', 'wt');
fprintf(fid, '%s\n', char(C));
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
