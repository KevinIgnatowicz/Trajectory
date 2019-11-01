function [t,x] = TraceXZ;
time= cputime;
variables;
%global Xin
%Trace 3D la trajectoire d'un objet sphérique dans un ecoulement
%Le bloc suivant lit le fichier CFD et créé les interpolants pour la
%vitesse

%global Points_U Points_V Points_W Xin
%{
A = textread(filename);

Pts = [A(:,18) A(:,19) A(:,20)]  ; % tableau des coordonnÃ©es des points du maillage

Mat_Vx =(A(:,2))./(A(:,1)); %vecteur colonne contenant les vitesses U en chaque point du maillage

Points_U = scatteredInterpolant(Pts, Mat_Vx, method); %interpolant de la vitesse U

Mat_Vy =(A(:,3))./(A(:,1)); %vecteur colonne contenant les vitesses V en chaque point du maillage

Points_V = scatteredInterpolant(Pts, Mat_Vy, method); %interpolant de la vitesse V

Mat_Vz =(A(:,4))./(A(:,1)); %vecteur colonne contenant les vitesses W en chaque point du maillage

Points_W = scatteredInterpolant(Pts, Mat_Vz, method); %interpolant de la vitesse W
%}
%La suite résoud le système différentiel et trace les graphes désirés
%vecteur du temps de simulation
xspan=[0.0:deltat:tmax] ;  

%ODE numerical solution:
 %options = odeset('OutputFcn',@odeplot);
[t, x]=ode23('Shimoi2', xspan, Xin);

%Mise des resultats sous forme graphique

taille=int16(tmax/deltat);
p1=x(:,1);
p2=x(:,3);
p3=x(:,5);

SH=zeros(taille,1);
for j = 1:taille
SH(j)=inShape(BWB,p1(j),p2(j),p3(j));
end
SH;
F=find(SH);
if Mode ==1
if size(F,1)==0
    
%Trajectoire seule
figure(1)
hold on 
plot3(x(:,1),x(:,3), x(:,5));
axis equal
grid on;
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
title('Trajectoires 3D');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
grid on;
hold off
end
end 

if Mode ==2
if size(F,1)==0

else
    n=F(1)-1;
    x_cut=p1(n);
    y_cut=p2(n);
    z_cut=p3(n);
    
figure(1)
hold on 
scatter3(x_cut,y_cut, z_cut, '.');
axis equal
title('Trajectoires 3D');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
grid on;
hold off
end
end

%resultats sous forme de tableaux


Z=zeros(taille,1);
Y=zeros(taille,1);
X=zeros(taille,1);
U=zeros(taille,1);
V=zeros(taille,1);
W=zeros(taille,1);
OMEGA=zeros(taille,1);
THETA=zeros(taille,1);
for i = 1:(taille)
    X(i,1)=x(i,1);
    Y(i,1)=x(i,3);
    Z(i,1)=x(i,5);
    U(i,1)=x(i,2);
    V(i,1)=x(i,4);
    W(i,1)=x(i,6);
    OMEGA(i,1)=x(i,8);
    THETA(i,1)=x(i,7);
end
X;
Y;
Z;
THETA;

%exporter les tableaux dans excel
filename = 'Magnus.xlsx';
xlswrite(filename,X,4,'V2')
xlswrite(filename,Z,4,'W2')
% xlswrite(filename,THETA,3,'P2')
% xlswrite(filename,U,3,'AG2')
% xlswrite(filename,W,3,'AI2')
% xlswrite(filename,OMEGA,3,'AK2')


%%%%%%%%%%%%%% Modif 2 bis : Sortie txt pour compatibilite Linux
%%% Precision a 10^-9
header1 = 'position-X';
header2 = 'position-Y';
fid=fopen('Magnus.txt','w');
fprintf(fid, [ header1 '\t' header2 '\n']);
fprintf(fid, '%11.9f \t %11.9f \n', [X Z]');
fclose(fid);

if size(F,1)==0 
    x_traj=X;
    y_traj=Y;
    z_traj=Z;
else
    x_traj=x_cut;
    y_traj=y_cut;
    z_traj=z_cut;
end
    Pr1= empreinte_trajectoire(x_traj,y_traj,z_traj, X_Plan);
    %Pr2= empreinte_trajectoire(x_traj,y_traj,z_traj, 2*X_Plan);

if size(Pr1,2)~=3
figure(2)
hold on
scatter(Pr1(1),Pr1(2),70,'.')
hold off
else
end
%{
if size(Pr2,2)~=3
figure(3)
hold on
scatter(Pr2(1),Pr2(2),70,'.')
hold off
else
end
%}


figure(3)
hold on
plot(-x(:,5),-x(:,6))
grid on
axis equal
hold off
%{
figure(4)
hold on
plot(x(:,1),x(:,5))
grid on
axis equal
hold off
%}

%e=cputime-time

%msg=msgbox(num2str(e),'Execution time')
%}
AA=[-x(:,5), -x(:,6)]
end

