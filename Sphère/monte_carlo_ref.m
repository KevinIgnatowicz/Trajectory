function [t,x] = monte_carlo_ref;
tic

%% Description

% (voir monte_carlo_main.m)

% Utilise la méthode Monte-Carlo pour calculer trajectoires moyennes et
% variances pour un grand nombre d'échantillons (800 000). L'objectif est
% d'utiliser ces valeurs comme référence pour la fonction monte-carlo_main.



%% Parametre de la simulation


%%% Physique

% Vitesse fluide
Vf=25;
% Angle de la plaque
angini=90;
% Temps de simulation
tmax=0.8;
% Erreur sur l'angle
deltang=2;


%%% Numerique

% Nombre d'échantillon
n=800000;
% Pas de temps
deltat=0.001;
% Nb d'itération en temps
taille=int16(tmax/deltat);



%% Initialisation des vecteurs [1]


% Vecteur du temps de simulation
xspan=[0.0:deltat:tmax]; 
T=zeros(taille,1);   % utilise par le solveur

% Vecteur parametre ideal
Xin=[0, 0, 0, 0, 90, 0, Vf]; % = [X0;u0;Z0;w0;tetha0;omega0;U_inf]

% Vecteurs solutions
Z=zeros(taille,1);
X=zeros(taille,1);
U=zeros(taille,1);
W=zeros(taille,1);
OMega=zeros(taille,1);
THETA=zeros(taille,1);

% Vecteurs moyenne
Zmoy=zeros(taille,1);
Xmoy=zeros(taille,1);
Umoy=zeros(taille,1);
Wmoy=zeros(taille,1);
THetamoy=zeros(taille,1);
OMegamoy=zeros(taille,1);

% Vecteurs min/max
Zmax=zeros(taille,1);
Xmax=zeros(taille,1);
Zmin=zeros(taille,1);
Xmin=zeros(taille,1);

% Variance
XVar = zeros(taille,1) ;
ZVar = zeros(taille,1) ;
% Vecteur de calcul intermediaire
Xtemp = zeros(taille,1) ;
Ztemp = zeros(taille,1) ;



%% Simulation sur n echantillons [1]

% - Generation Monte-Carlo          OK
% - Resolution numerique (ode23)    OK
% - Vecteurs moyennes               OK
% - Vecteurs variances              OK
% - Trajectoires extremales         OK


%%% Boucle de calcul
for j=1:n    
    j % affiche l'iteration en cours
    
    %%% Generation du vecteur inital
    
    % Signe aleatoire
    sgn1 = -1 + 2*randi([0 1],1);
    sgn2 = -1 + 2*randi([0 1],1);
    sgn3 = -1 + 2*randi([0 1],1);
    
    theta0 = angini - deltang/2 + (deltang*rand(1));    % erreur : +/- deltang
    omega0=rand(1)*0.5*sgn1;                            % erreur : 0 à 1 rad/s
    u0=rand(1)*(3*Vf/100)*sgn2;                         % erreur : +/-0.03.Vf
    w0=rand(1)*(3*Vf/100)*sgn3;                         % erreur : +/-0.03.Vf
    
    % Enregistrement de Theta0_echantillon
    Plage(j,1)=theta0;

    Xin=[0, u0, 0, w0, theta0, omega0,Vf]; % = [X0;u0;Z0;w0;tetha0;omega0]


    %%% Resolution numerique
    [t, x]=ode45('Shimoi2', xspan, Xin);
    % Remplacer la ligne du dessus par le bloc suivant pour tracer le
    % vecteur solution.
    %figure(1);
    %options = odeset('OutputFcn',@odeplot); 
    %[t, x]=ode45('motion_plate_Shimoi', xspan, Xin,options);
    
    X(1:taille,1) = x(1:taille,1);
    Z(1:taille,1) = x(1:taille,3);
    U(1:taille,1) = x(1:taille,2);
    W(1:taille,1) = x(1:taille,4);
    OMega(1:taille,1) = x(1:taille,6);
    THETA(1:taille,1) = x(1:taille,5);
    T=t;

    
    %%% Calcul des moyennes
    Xmoy=(Xmoy+X);
    Zmoy=(Zmoy+Z);
    Umoy=(Umoy+U);
    Wmoy=(Wmoy+W);
    THetamoy=(THetamoy+THETA);
    OMegamoy=(OMegamoy+OMega);
    
    
    %%% Variance Empirique
    Xtemp = Xtemp + X.^2 ;
    Ztemp = Ztemp + Z.^2 ;
    XVar = Xtemp/j - (Xmoy/j).^2 ;
    ZVar = Ztemp/j - (Zmoy/j).^2 ;
    
    %%% Calcul des Trajectoires extremales
    for k=1:taille
        if j==1
            Zmax(k,1)=Z(k,1);
            Xmax(k,1)=X(k,1);
            Zmin(k,1)=Z(k,1);
            Xmin(k,1)=X(k,1);    
        else
            Zmax(k,1)=max(Zmax(k,1),Z(k,1));
            Zmin(k,1)=min(Zmin(k,1),Z(k,1));
            Xmax(k,1)=max(Xmax(k,1),X(k,1));
            Xmin(k,1)=min(Xmin(k,1),X(k,1));
        end
    end
    
    
end



%% Sortie Donnes

% - Export txt                          a completer
% - Export excel                        a completer


%%% Sortie texte, compatible Linux
% Precision a 10^-9

fid=fopen('Monte_Carlo_800000_X.txt','w');
fprintf(fid, '%11.9f \n', ...
              [Xmoy/n]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_Z.txt','w');
fprintf(fid, '%11.9f \n', ...
              [Zmoy/n]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_VX.txt','w');
fprintf(fid, '%11.9f \n', ...
              [XVar]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_VZ.txt','w');
fprintf(fid, '%11.9f \n', ...
              [ZVar]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_XMIN.txt','w');
fprintf(fid, '%11.9f \n', ...
              [Xmin]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_ZMAX.txt','w');
fprintf(fid, '%11.9f \n', ...
              [Xmax]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_ZMIN.txt','w');
fprintf(fid, '%11.9f \n', ...
              [Zmin]');
fclose(fid);

fid=fopen('Monte_Carlo_800000_ZMAX.txt','w');
fprintf(fid, '%11.9f \n', ...
              [Zmax]');
fclose(fid);



toc


end