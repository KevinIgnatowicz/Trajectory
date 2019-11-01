function [t,x] = monte_carlo_ex;
tic

%% Description

% (voir monte_carlo_main.m)

% Illustration de la Méthode de Monte-Carlo :

% 100 Trajectoires sont générés, en générant aléatoirement les paramètres initiaux
% dans le spectre d'incertitude donné. On calcule alors la moyenne et variance
% empirique pour chaque pas de temps, pour simuler une expérience de
% laboratoire en prenant en compte les incertitudes de mesures.

% Toutes les trajectoires obtenues par les échantillons générés
% aléatoirement sont tracées, ainsi que la moyenne, les trajectoires extrémales et
% l'écart type à chaque pas de temps en Z.
% Le pas de temps est choisit volontairement élevé pour ne pas surcharger
% les graphes.
% Les axes sont uniformisés pour une meilleure représentation : Ils
% peuvent être à revoir selon des cas particuliers



%% Parametre de la simulation


%%% Physique

% Vitesse fluide
Vf=25;
% Angle de la plaque
angini=0;
% Temps de simulation
tmax=0.8;
% Erreur sur l'angle
deltang=180;


%%% Numerique

% Nombre d'échantillon
n=100;
% Pas de temps
deltat=0.02;
% Nb d'iteration en temps
taille=int16(tmax/deltat);



%% Initialisation des vecteurs


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
% - Trajectoires extremales         OK
% - Export graphique                OK


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

    Xin=[0, u0, 0, w0, theta0, omega0,Vf]; % = [X0;u0;Z0;w0;tetha0;omega0;U_inf]


    %%% Resolution numerique
    [t, x]=ode45('Shimoi2', xspan, Xin);
    
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

  
  %%% Export Graphique
  figure(1);
  axis('equal')
  axis([0 max(Xmin) -3.5 3])
  plot(X(:,1),Z(:,1));
  title('Path examples');
  xlabel('X Path (meter)')
  ylabel('Z Path (meter)')
  axis('equal')
  axis([0 max(Xmin) -3.5 3])
  grid on;
  hold on;


    
    
end



%% Enveloppe de trajectoire

% Interpolation des courbes min/max pour tracer l'enveloppe des courbes
% extrémales. Pas la meilleure méthode (limité par matlab) mais produit
% des résultats graphiques exploitables.

% Les courbes extrémales ne sont pas tout à fait exactes. Un encadrement
% plus précis nécéssiterais une interpolations de chaque trajectoire sur la même grille.
% Il est sans doute préférable de trouver une meilleure fonction pour créer
% l'enveloppe.

% grille d'interpolation
xq= 0:0.008:max(Xmin) ; 
zq= -3:0.003:3 ; 

% courbe interpolée
v1q = interp1(Xmin,Zmax,xq);
v2q = interp1(Xmin,Zmin,xq);

%#create continuous x value array for plotting
A=[xq,fliplr(xq)];                
B=[v1q,fliplr(v2q)];


%% Sortie Graphique

% - Figure 1 : Trajectoires             OK
% - Figure 2 : Moyenne et Enveloppe     OK
%  -Figure 2 : Moyenne et écart type    OK


% Trajectoire moyenne et enveloppe
figure(2);
  h = fill(A,B,'b');
  set(h,'facealpha',.1)
  hold on;
  plot(Xmoy(:,1)/n,Zmoy(:,1)/n,'color','r');
  title('Mean and extremal paths');
  xlabel('X Path (meter)')
  ylabel('Z Path (meter)')
  axis('equal')
  axis([0 max(Xmin) -3.5 3])
  grid on;
  
figure(3);
  errorbar(Xmoy(:,1)/n,Zmoy(:,1)/n,sqrt(XVar));
  hold on;
  plot(Xmoy(:,1)/n,Zmoy(:,1)/n,'color','r');
  title('Mean path and X standard deviation');
  xlabel('X Path (meter)')
  ylabel('Z Path (meter)')
  axis('equal')
  axis([0 max(Xmin) -3.5 3])
  grid on;
 
 toc


end