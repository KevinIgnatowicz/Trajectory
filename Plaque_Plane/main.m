function [t,x ] = main; 
time=cputime;
variables_plaque;
%Le bloc suivant lit le fichier CFD et créé les interpolants pour la
%vitesse

global Points_U Points_V Points_W  Xin


A = csvread(filename,1,0);

Pts = [A(:,24) A(:, 25) A(:, 26)]  ; % tableau des coordonnees des points du maillage

Mat_Vx =(A(:,2))./(A(:,1)); %vecteur colonne contenant les vitesses U en chaque point du maillage

Points_U = scatteredInterpolant(Pts, Mat_Vx, method); %interpolant de la vitesse U

Mat_Vy =(A(:,3))./(A(:,1)); %vecteur colonne contenant les vitesses V en chaque point du maillage

Points_V = scatteredInterpolant(Pts, Mat_Vy, method); %interpolant de la vitesse V

Mat_Vz =(A(:,4))./(A(:,1)); %vecteur colonne contenant les vitesses W en chaque point du maillage

Points_W = scatteredInterpolant(Pts, Mat_Vz, method); %interpolant de la vitesse W

%La suite résoud le système différentiel et trace les graphes désirés
%this function solve ODE for rotation angle using quaternion
 
%note le pas de temps est variable dans ode23

%vecteur du temps de simulation
tspan=[0.0:deltat:tmax]; 

%ODE numerical solution:
 %options = odeset('OutputFcn',@odeplot);
  

[t, x]=ode23('Plaque_plane_6DDL', tspan, Xin);


%tracé de la trajectoire
taille=size(x,1);

p1=x(:,8);
p2=x(:,10);
p3=x(:,12);

SH=zeros(taille,1);
for j = 1:taille
SH(j)=inShape(BWB,p1(j),p2(j),p3(j));
end
SH;
F=find(SH);
if size(F,1)==0 
    
%Trajectoire seule

figure(1)
hold on 
plot3(x(:,8),x(:,10), x(:,12));
axis equal
grid on;
saveas(gcf,'Trajectoire_NOV18_Richards','fig');
hold off

else
    n=F(1)-1;
    x_cut=p1(1:n);
    y_cut=p2(1:n);
    z_cut=p3(1:n);
    
figure(1)
hold on 
plot3(x_cut,y_cut, z_cut);
axis equal
grid on;
saveas(gcf,'Trajectoire_NOV18_Richards','fig');
hold off
end

% les 2 graphes suivants permettent au choix de tracer d'autres paramètres
%{
figure(3)
hold on
plot(t,x(:,10));
title('top');
xlabel('X')
ylabel('Y')
grid on;
hold off

figure(4)
hold on
plot(t,x(:,12));
title('side');
xlabel('X')
ylabel('Z')
grid on;
hold off
%}

%resultats sous forme de tableaux

Z=zeros(taille,1);
Y=zeros(taille,1);
X=zeros(taille,1);
U=zeros(taille,1);
V=zeros(taille,1);
W=zeros(taille,1);
P=zeros(taille,1);
Q=zeros(taille,1);
R=zeros(taille,1);
Phi=zeros(taille,1);
Theta=zeros(taille,1);
Psi=zeros(taille,1);
for i = 1:(taille)
    X(i,1)=x(i,8);
    Y(i,1)=x(i,10);
    Z(i,1)=x(i,12);
    U(i,1)=x(i,9);
    V(i,1)=x(i,11);
    W(i,1)=x(i,13);
    P(i,1)=x(i,1);
    Q(i,1)=x(i,2);
    R(i,1)=x(i,3);
    Phi(i,1)=x(i,14);
    Theta(i,1)=x(i,15);
    Psi(i,1)=x(i,16);
end
X;
Y;
Z;
U;
V;
W;
P;
Q;
R;
Phi;
Theta;
Psi;
%{
Vrx=zeros(taille,1);
Vry=zeros(taille,1);
Vrz=zeros(taille,1);
Quat=zeros(taille,4);
UB=zeros(taille,3);


for i=1:taille
    Vrx(i,1)= velocity(x(i,8),x(i,10),x(i,12),1)-U(i,1);
    Vry(i,1)= velocity(x(i,8),x(i,10),x(i,12),2)-V(i,1);
    Vrz(i,1)= velocity(x(i,8),x(i,10),x(i,12),3)-W(i,1);
    [Quat(i,1),Quat(i,2),Quat(i,3),Quat(i,4)]=Euler2quat(Phi(i,1),-Theta(i,1),-Psi(i,1));
    VB=Rq(Quat(i,1),Quat(i,2),Quat(i,3),Quat(i,4))*[Vrx(i,1);Vry(i,1);Vrz(i,1)];
    UB(i,1)=VB(1);
    UB(i,2)=VB(2);
    UB(i,3)=VB(3);
 end
    

Md_X=zeros(taille,1);
Md_Y=zeros(taille,1);
Md_Z=zeros(taille,1);
if choix_moment == 1
    for i=1:(taille)
        Md_X(i,1)=-0.01*0.5*rhoair*(sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)+sqrt(P(i,1).^2+Q(i,1).^2+R(i,1).^2)*l/2)*L*l.^3.*P(i,1);
        Md_Y(i,1)=-0.185*0.5*rhoair*(sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)+sqrt(P(i,1).^2+Q(i,1).^2+R(i,1).^2)*l/2)*L*l.^3.*Q(i,1);
        Md_Z(i,1)=-0.2*0.5*rhoair*(sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)+sqrt(P(i,1).^2+Q(i,1).^2+R(i,1).^2)*l/2)*L*l.^3.*R(i,1);
    end
end
if choix_moment == 2
    W0=zeros(taille,1);
    W1=zeros(taille,1);
    for i=1:(taille)
        if L/l == 1
W0(i,1)=0.64*sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)/L; %W0 = 0.32*2V/L plaque carrée
W1(i,1)=0.64*sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)/L;
end
if L/l == 2
W0(i,1)=(-1.25*(e/l)+0.5163)*2*sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)/l; %W0 = 0.48*2V/l plaque rectangulaire
W1(i,1)=0.8*sqrt(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)/L;
end
   Md_X(i,1)=0.5*rhoair*(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)*l*L*l.*sign(P(i,1)).*Tachikawa(abs(P(i,1))./W0(i,1),L,l);
   Md_Y(i,1)=0.5*rhoair*(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)*L*L*l.*sign(Q(i,1)).*Tachikawa(abs(Q(i,1))./W1(i,1),L,l);
   Md_Z(i,1)=0.5*rhoair*(UB(i,1).^2+UB(i,2).^2+UB(i,3).^2)*L*L*e.*sign(R(i,1)).*Tachikawa(abs(R(i,1))./W0(i,1),L,l);
    end
end
 Md_X;
 Md_Y;
 Md_Z;
%} 
%{
figure(5)
hold on
plot(t,x(:,14))
title('Phi');
xlabel('t (s)')
ylabel('Phi')
grid on;
hold off

figure(6)
hold on
plot(t,x(:,15));
title('Theta');
xlabel('t (s)')
ylabel('Theta')
grid on;
hold off

figure(7)
hold on
plot(t,x(:,16));
title('Psi');
xlabel('t (s)')
ylabel('Psi')
grid on;
hold off

figure(8)
hold on
plot(t,x(:,2));
title('Q');
xlabel('t (s)')
ylabel('Q')
grid on;
hold off

figure(2)
hold on
plot(x(:,8),x(:,2));
title('Q');
xlabel('X (m)')
ylabel('Q')
grid on;
hold off

figure(3)
hold on
plot(x(:,8),Md_Y);
title('Md_Y');
xlabel('X (m)')
ylabel('Md_Y')
grid on;
hold off
%}


if size(F,1)==0 
    x_traj=X;
    y_traj=Y;
    z_traj=Z;

else
    x_traj=x_cut;
    y_traj=y_cut;
    z_traj=z_cut;
end
global Pr
    Pr= empreinte_trajectoire(x_traj,y_traj,z_traj, X_Plan);

if Pr(1)~='N'
figure(2)
hold on
scatter(Pr(1),Pr(2),30,'s','filled')
saveas(gcf,'Empreintes_NOV18_Richards','fig');
hold off
else
end

%timef=cputime-time

%msg=msgbox(num2str(timef),'Execution time')

end

