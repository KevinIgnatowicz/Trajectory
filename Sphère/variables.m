filename = 'ONERAM6_M05.txt'; %Fichier CFD donnant l'écoulement
B = csvread('surface_flow.csv',1,0);
BWB=alphaShape(B(:,2), B(:,3), B(:,4),2);

method = 'linear'; % Méthode d'interpolation pour les vitesses
Mode=1 % 1 : trajectoires, 2 : collection gouttelettes

X_Plan = 20; %Abscisse (m) du plan de normale x dans lequel on veut visualiser l'intersection avec les trajectoires

Xin=[0,0,0,0,0,-1,0,0]; %[X0;u0;Y0;v0;Z0;w0;tetha0;omega0];

%temps de simulation
tmax=0.3;

%pas de temps
deltat=0.0001;

g=9.81; %gravité

rhoShimoi=1000;%57.25*0.45359237/(0.3048^3); % 917.0570 masse vol de la glace kg/m^3 annoncee par Koji Shimoi
rhoPlastique=1120; % masse volumique du plastique(kg/m^3)


%%constante pour une sphere
dsp=0.000100 ;% diametre


muair=1.8e-5; % viscosite de l'air
rhoair= 1.2237; %Rho2D(filename, x(1),x(3),method); % densite de l'air (kg/m^3)



