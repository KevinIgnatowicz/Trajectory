function [t,x] = monte_carlo_main;
tic

%% Description

% Utilise la méthode Monte-Carlo pour simuler les incertitudes
% expérimentales sur le code de calcul de plaque plane dans un
% écoulement uniforme.

% L'équation différentielle est résolue numériquement par la fonction
% 'Shimoi2' (précédement 'motion_plate_Shimoi') 

% Le principe de la méthode est de générer un nombre important
% d'échantillons de paramètres initiaux, assortis d'une erreur
% aléatoire comprise dans l'amplitude des incertitudes (définie par l'utilisateur).

% Un critère de précision epsilon est définit en première partie du code. 
% La fonction effectue 100 premières simulations et utilise les résultats
% obtenus pour prévoir le nombre d'itérations nécessaire pour satisfaire
% les relations suivantes :
% |g(n+1) - g(n)| < epsilon
% |s(n+1) - s(n)| < epsilon
% Avec g la trajectoire moyenne et s la variance empirique calculé après n
% simulations

% Si les résultats de références ont été calculés (exécution de la fonction
% 'monte_carlo_ref' avec les mêmes paramètres), la distance à la solution
% de référence est également calculée (courbes d'erreur)

% Les différents résultats obtenus sont exportés sous format txt
 


% fonction modifiee depuis la version originale : 'monthe_carlo1'
% Principaux ajouts :
% - Correctiond de la methode (variable aleatoire decorellées)
% - Ajouts des calculs de variation et majorations
% - Prévision du nombre d'itérations nécessaires pour un critère d'erreur donné
% - Sorties graphiques et texte supplementaires
% - Modification du solveur pour prendre Vf comme variable d'entrée


% Description


%% Parametre de la simulation


%%% Physique

% Vitesse fluide
Vf=25;
% Angle de la plaque
angini=0;
% Temps de simulation
tmax=0.8;
% Erreur sur l'angle
deltang=4;


%%% Numerique

% Nombre d'�chantillon
n=101;
% Pas de temps
deltat=0.001;
% Nb d'iteration en temps
taille=int16(tmax/deltat);
% Critere d'erreur max.
epsilon = 0.1;
% Compare les résultats obtenus avec les valeurs de références si les
% fichiers sont trouvés
comp = 0;

comp = comp + exist('Monte_Carlo_800000_X.txt','file') ...
            + exist('Monte_Carlo_800000_Z.txt','file') ...
            + exist('Monte_Carlo_800000_VX.txt','file') ...
            + exist('Monte_Carlo_800000_VZ.txt','file') ;



%% Initialisation des vecteurs [1]


% Vecteur du temps de simulation
xspan=[0.0:deltat:tmax]; 
T=zeros(taille,1);   % utilise par le solveur

% Vecteur parametre ideal
Xin=[0, 0, 0, 0, 90, 0, Vf]; % = [X0;u0;Z0;w0;tetha0;omega0;U_infini]

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

% Vecteurs min/max (a utiliser)
Zmax=zeros(taille,1);
Xmax=zeros(taille,1);
Zmin=zeros(taille,1);
Xmin=zeros(taille,1);

% Vecteurs facteur de majoration 
F=zeros(20,1);  % Moyenne en X
G=zeros(20,1);  % Moyenne en Z

L=zeros(20,1);  % Variance en X
M=zeros(20,1);  % Variance en Z



%% Initialisation des vecteurs [2]


%%% Variation de la moyenne

% Moyenne a l'iteration precedente
XOLDmoy=zeros(taille,1);
ZOLDmoy=zeros(taille,1);

% Variation de la moyenne en norme infini
XVmoy=zeros(n,1);
ZVmoy=zeros(n,1);

% Majoration de la variation de moyenne
XVmaj=zeros(n,1);
ZVmaj=zeros(n,1);


%%% Vecteur Erreur relative

% Solution ideale
XT=zeros(taille,1);
ZT=zeros(taille,1);
UT=zeros(taille,1);
WT=zeros(taille,1);
THetaT=zeros(taille,1);
OMegaT=zeros(taille,1);

% Distance |Moy.echantillon - Sol.ideale| en norme inifni
XTmoy=zeros(n,1);
ZTmoy=zeros(n,1);


%%% Vecteur variance

% Variance
XVar = zeros(taille,1) ;
ZVar = zeros(taille,1) ;


% Vecteur de calcul intermediaire (utilise pour le calcul de la variance)
Xtemp = zeros(taille,1) ;
Ztemp = zeros(taille,1) ;


%%% Vecteur Erreur relative

% Variance Solution ideale (800 000 iterations)
XTVar=zeros(taille,1);
ZTVar=zeros(taille,1);


% Distance |VAR - TVAR| en norme inifni
XDVar=zeros(n,1);
ZDVar=zeros(n,1);

% Variance a l'iteration precedente
XOLDVar=zeros(taille,1);
ZOLDVar=zeros(taille,1);

% Variation de la variance en norme infini
XVVar=zeros(n,1);
ZVVar=zeros(n,1);




%% Lecture des trajectoires ideales

% Lecture des moyenne et variances à 800 000 simulations si les fichiers existent

if comp == 8
    disp('Fichiers trouvés. Lecture des valeurs de référence')
    
    fileID = fopen('Monte_Carlo_800000_X.txt','r') ;
    formatSpec = '%f';
    sizeXT = [1 Inf];
    XT0 = fscanf(fileID,formatSpec,sizeXT) ;
    fclose(fileID);

    fileID = fopen('Monte_Carlo_800000_Z.txt','r') ;
    formatSpec = '%f';
    sizeZT = [1 Inf];
    ZT0 = fscanf(fileID,formatSpec,sizeZT) ;
    fclose(fileID);

    XT = XT0.' ;
    ZT = ZT0.' ;

    % Variance

    fileID = fopen('Monte_Carlo_800000_VX.txt','r') ;
    formatSpec = '%f';
    sizeXT = [1 Inf];
    XT0 = fscanf(fileID,formatSpec,sizeXT) ;
    fclose(fileID);

    fileID = fopen('Monte_Carlo_800000_VZ.txt','r') ;
    formatSpec = '%f';
    sizeZT = [1 Inf];
    ZT0 = fscanf(fileID,formatSpec,sizeZT) ;
    fclose(fileID);

    XTVar = XT0.' ;
    ZTVar = ZT0.' ;
end

j=0;
flag = false ;


%% Simulation sur n echantillons [1]

% - Generation Monte-Carlo          OK
% - Resolution numerique (ode23)    OK
% - Vecteurs moyennes               OK
% - Trajectoires extremales         non utilisées ici


%%% Boucle de calcul

% Le nombre réel d'échantillons requis n'est connue qu'à la centième itération.
% Ce nombre est calculé en fonction du critère d'arret et par curve-fitting sur 20 valeurs.

while flag == false
    j = j+1 % affiche l'iteration en cours
    n       % affiche l'iteration max. prevue

    %%% Generation du vecteur d'erreur
    
    % Signe aleatoire
    sgn1 = -1 + 2*randi([0 1],1);
    sgn2 = -1 + 2*randi([0 1],1);
    sgn3 = -1 + 2*randi([0 1],1);
    
    theta0 = angini - deltang/2 + (deltang*rand(1));    % erreur : +/- deltang/2
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


    
%% Simulation sur n echantillons [2]

% - Variation de la moyenne             OK
% - Variance Empirique                  OK
% - Erreur /sol. ideale                 OK
% - Majoration des variations           OK



    %%% XVmoy  : Variation de la moyenne
    if j==1
        XVmoy(j,1)=0;
        ZVmoy(j,1)=0;  
    else      
        XVmoy(j,1)=max( abs(Xmoy/j - XOLDmoy/(j-1)) );
        ZVmoy(j,1)=max( abs(Zmoy/j - ZOLDmoy/(j-1)) );       
    end
    XOLDmoy = Xmoy;
    ZOLDmoy = Zmoy;


    
    %%% Variance Empirique
    Xtemp = Xtemp + X.^2 ;
    Ztemp = Ztemp + Z.^2 ;
    XVar = Xtemp/(j^2) - (Xmoy/j).^2 ;
    ZVar = Ztemp/(j^2) - (Zmoy/j).^2 ;
    XVarMax(j,1) = max( abs(XVar) );
    ZVarMax(j,1) = max( abs(ZVar) );
    
    
    %%% XVVar  : Variation de la variance
    if j==1
        XVVar(j,1)=0;
        ZVVar(j,1)=0;
    else      
        XVVar(j,1)=max( abs(XVar - XOLDVar) );
        ZVVar(j,1)=max( abs(ZVar - ZOLDVar) );
    end
    XOLDVar = XVar;
    ZOLDVar = ZVar;
    
    %%% Les calculs suivant ne sont effectuées que si les fichiers de références
    %%% sont détectées
  
    if comp == 8

        %%% XTmoy : Distance entre moyenne et sol.ideale
        XTmoy(j,1) = max( abs( XT - Xmoy/j ) ) ;
        ZTmoy(j,1) = max( abs( ZT - Zmoy/j ) ) ;
   
        %%% XDVar : Distance entre variance et sol.ideale
        XDVar(j,1) = max( abs( XTVar - XVar ) ) ;
        ZDVar(j,1) = max( abs( ZTVar - ZVar ) ) ;

    end

    
    %%% Calcul du coefficient d'erreur sur les itérations 80..100.
    
    % Le facteur d'erreur est calculé pour les moyennes et variances en X
    % et en Z, donnant 4 valeurs du nombre d'itérations nécessaires
    % pour le même critère epsilon.
    % Le maximum des 4 valeurs est ensuite choisi.

    
    if j==80
        F(1,1) = XVmoy(j,1)*j ;
        G(1,1) = ZVmoy(j,1)*j ;    
        L(1,1) = XVVar(j,1)*j ;
        M(1,1) = ZVVar(j,1)*j ;
    elseif j>80 && j<100
        F(j-79,1) = XVmoy(j,1)*j ;
        G(j-79,1) = ZVmoy(j,1)*j ;
        L(j-79,1) = XVVar(j,1)*j ;
        M(j-79,1) = ZVVar(j,1)*j ;
    elseif j==100
        F(j-79,1) = XVmoy(j,1)*j ;
        G(j-79,1) = ZVmoy(j,1)*j ;
        L(j-79,1) = XVVar(j,1)*j ;
        M(j-79,1) = ZVVar(j,1)*j ;
        
        fmaj = max(F) ;
        gmaj = max(G) ;
        lmaj = max(L) ;
        mmaj = max(M) ;
        
        n1 = ceil(fmaj/epsilon) ;
        n2 = ceil(gmaj/epsilon) ;
        n3 = ceil(lmaj/epsilon) ;
        n4 = ceil(mmaj/epsilon) ; 
        n  = max([n1,n2,n3,n4]) ;
    end
    
    % Arrêt des itérations si le nombre requis est atteint
    if j>n-1 && j>99
       flag = true;
    end


    
end





%% Post-traitement

% - Comparaison avec majoration theorique               OK

%%% Majoration theorique des variations de la moyenne
for k=1:n
    XVmaj(k,1) = fmaj/k ;
    ZVmaj(k,1) = gmaj/k ;
end

%%% Majoration theorique des variations de la variance
for k=1:n
    XVVarmaj(k,1) = lmaj/k ;
    ZVVarmaj(k,1) = mmaj/k;
end



%% Sortie graphique

% Variation de la moyenne (X) avec majoration / echelle log
figure(1);
  plot(log(XVmoy(:,1)),'color',[0.5,0.5,0.5]);
  hold on;
  plot(log(XVmaj(:,1)),'color','k')
  title('V_f = 25, Theta_0 = 0');
  xlabel('Nombre de simulations')
  ylabel('Variation de moyenne (X)')
  grid on;

  
% Variation de la moyenne (Z) avec majoration / echelle log
figure(2);
  plot(log(ZVmoy(:,1)),'color',[0.5,0.5,0.5]);
  hold on;
  plot(log(ZVmaj(:,1)),'color','k')
  title('V_f = 25, Theta_0 = 0');
  xlabel('Nombre de simulations')
  ylabel('Variation de moyenne (Z)')
  grid on;

  
% Variation de variance (X)
figure(7);
  plot(log(XVVarmaj(:,1)),'color','k');
  hold on;
  plot(log(XVVar(:,1)),'color',[0.5,0.5,0.5])
  title('V_f = 25, Theta_0 = 0');
  xlabel('Nombre de simulations')
  ylabel('Variation de variance (X)')
  grid on;
  
  
% Variation de variance (Z)
figure(8);
  plot(log(ZVVarmaj(:,1)),'color','k');
  hold on;
  plot(log(ZVVar(:,1)),'color',[0.5,0.5,0.5])
  title('V_f = 25, Theta_0 = 0');
  xlabel('Nombre de simulations')
  ylabel('Variation de variance (Z)')
  grid on;

  
%%% Les courbes suivantes ne sont tracées que si les fichiers de références
%%% sont détectées
  
if comp == 8
    
% Erreur en norme infinie (X)
figure(3);
  plot(log(XTmoy(:,1)),'color','k');
  title('X-Path error, \Delta\theta = 2');
  xlabel('Number of simulations')
  ylabel('Uniform norm error (meter)')
  grid on;
 
  
% Erreur en norme infinie (Z)
figure(4);
  plot(log(ZTmoy(:,1)),'color','k');
  title('Z-Path error, \Delta\theta = 2');
  xlabel('Number of simulations')
  ylabel('Uniform norm error (meter)')
  grid on;
  
  
% Erreur en norme infinie (VX)
figure(5);
  plot(log(XDVar(:,1)),'color','k');
  title('X-Variance error, \Delta\theta = 2');
  xlabel('Number of simulations')
  ylabel('Uniform norm error (meter)')
  grid on;
 
  
% Erreur relative en norme infinie (VZ)
figure(6);
  plot(log(ZDVar(:,1)),'color','k');
  title('Z-Variance error, \Delta\theta = 2');
  xlabel('Number of simulations')
  ylabel('Uniform norm error (meter)')
  grid on;
  
end




  
%% Sortie Donnes

% - Trajectoires moyennes et variances      OK
% - Variation de moyenne et majoration      OK
% - Variation de variance et majoration     OK
% - Erreur                                  OK

%%% Sortie texte, compatible Linux
% Precision a 10^-9


% Moyenne et variance
header1 = 'position X';
header2 = 'vitesse  U';
header3 = 'position Z';
header4 = 'position W';
header5 = 'Angle thet';
header6 = 'Angle alph';
header7 = 'Variance X';
header8 = 'Variance Z';

fid=fopen('Monte_Carlo_Moyenne.txt','w');
fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t'...
               header4 '\t' header5 '\t' header6 '\t'...
               header7 '\t' header8 '\n']);
fprintf(fid, '%11.9f \t %11.9f \t %11.9f \t %11.9f \t %11.9f \t %11.9f \t %11.9f \t %11.9f \n', ...
              [Xmoy/n Umoy/n Zmoy/n Wmoy/n THetamoy/n OMegamoy/n XVar ZVar]');
fclose(fid);


% Variation de moyenne et majoration
header1 = 'VarMoy X';
header2 = 'VarMajMoy X';
header3 = 'VarMoy Z';
header4 = 'VarMajMoy Z';

fid=fopen('Monte_Carlo_VariationMoy.txt','w');
fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\n']);
fprintf(fid, '%11.9f \t %11.9f \t %11.9f \t %11.9f \n', ...
              [XVmoy XVmaj ZVmoy ZVmaj]');
fclose(fid);


% Variation de variance et majoration
header1 = 'Var X';
header2 = 'VarVar X';
header3 = 'Var Z';
header4 = 'VarVar Z';

fid=fopen('Monte_Carlo_VariationVar.txt','w');
fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\n']);
fprintf(fid, '%11.9f \t %11.9f \t %11.9f \t %11.9f \n', ...
              [XVVar XVVarmaj ZVVar ZVVarmaj]');
fclose(fid);


% Erreur
header1 = 'Erreur Xmoy';
header2 = 'Erreur Zmoy';
header3 = 'Erreur XV';
header4 = 'Erreur ZV';

fid=fopen('Monte_Carlo_Erreur.txt','w');
fprintf(fid, [ header1 '\t' header2 '\t' header3 '\t' header4 '\n']);
fprintf(fid, '%11.9f \t %11.9f \t %11.9f \t %11.9f \n', ...
              [XTmoy ZTmoy XDVar ZDVar]');
fclose(fid);

toc



end