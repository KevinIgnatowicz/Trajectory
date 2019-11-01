filename='Clip_BWB08.csv'; %Fichier CFD contenant les données de l'écoulement
B = csvread('Surface_BWB.csv',1,0); %Fichier contenant les coordonnées de la surface de l'avion
BWB=alphaShape(B(:,24), B(:,25), B(:,26),2); %Définition de l'AlphaShape 
method = 'linear'; %linear, nearest ou natural. Méthode pour l'interpolation des vitesses
choix_moment = 1 ; % 1 : Richards; 2 : Tachikawa
X_Plan = 20; %Abscisse (m) du plan de normale x dans lequel on veut visualiser l'intersection avec les trajectoires
SR=1; %1,2 Surface Ratio de la plaque
L=0.1;  % Longueur de la plaque (m)   
l=L/SR ; %Largeur de la plaque (m), calculée automatiquement avec le SR    
e=0.002; %épaisseur de la plaque (m)
g=9.81;  %gravité (m/s2)
rhoShimoi=917; %masse volumique de la glace (kg/m3)
mass=rhoShimoi*L*l*e; %masse de la plaque (kg)
rhoair=1.195; %masse volumique de l'air (kg/m3)

%angles = [30;20;30]; %angles d'Euler initiaux (repère avion, rotations autour de z, y et x) en degrés
%[q10,q20,q30,q40]=Euler2quat(angles(1)*pi/180,-angles(2)*pi/180,-angles(3)*pi/180); %quaternions initiaux
%Conditions initiales : 
 %  [P0 Q0 R0 q10 q20 q30 q40  x0 u0  y0 v0 z0 w0      phi0               theta0             psi0]
%Xin=[ 0, 0, 0,q10,q20,q30,q40, 0, 0,  6, 0, 1, 0, angles(1)*pi/180, -angles(2)*pi/180, -angles(3)*pi/180];

%temps de simulation
tmax=0.55;
%nous forçons l'écriture à un intervalle deltat
%pas de temps
deltat=0.001;

%calcul de la matrice d'inertie de la plaque
  
Ixx=1/12*mass*(l^2+e^2);
Iyy=1/12*mass*(L^2+e^2);
Izz=1/12*mass*(l^2+L^2);
Iplate=diag([Ixx; Iyy; Izz]); % matrice d'inertie de la plaque
invIplate=inv(Iplate);  % inverse de la matrice d'inertie

% Choix de la base de donnée de Richards pour avoir le CN correspondant au
% SR choisi
if SR == 1
    data_CN = importdata('sr1.txt');
elseif SR==2
    data_CN = importdata('sr2.txt');
elseif SR == 4
    data_CN = importdata('sr4.txt');
end

alp = data_CN(:,1);
bet = data_CN(:,3);

Cn = data_CN(:,2);

Vq = scatteredInterpolant(alp,bet, Cn); %interpolant du CN