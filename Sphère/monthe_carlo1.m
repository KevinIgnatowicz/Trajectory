function [t,x] = monthe_carlo1;
tic
%%FLAT PLATE MOTION%%

%Trace dans le plan XZ la trajectorie d'une plaque plane dans un �coulement
%uniforne de vitesse Vf.

Vf=25;

%nombre d'�chantillon
n=3000;
%temps de simulation
tmax=0.8;
%pas de temps
deltat=0.001;

%initialisation des vecteurs.
taille=int16(tmax/deltat);
Zmoy=zeros(taille,1);
Xmoy=zeros(taille,1);

Umoy=zeros(taille,1);
Wmoy=zeros(taille,1);

THetamoy=zeros(taille,1);
OMegamoy=zeros(taille,1);

Zmax=zeros(taille,1);
Xmax=zeros(taille,1);

Zmin=zeros(taille,1);
Xmin=zeros(taille,1);

Repartition=zeros(taille,1);
Plage=zeros(n,1);
it=0;



for j=1:n

    it=j;
%signe al�atoire
sgne=rand(1) ;
if sgne>0.5;
    sgn1=-1;
else
    sgn1=1;
end
    
     
%%%%angle
angini=90;
%%%%%

%Delatangle
deltang=4;
%theta0
alea=rand(1);
theta0=angini-deltang/2+(deltang*alea);

Plage(j,1)=theta0;


%vitesse angulaire initiale
omega0=rand(1)*0.5*sgn1;
%vitesse en x initiale
u0=rand(1)*(3*Vf/100)*sgn1;
%vitesse en z initiale
w0=rand(1)*(3*Vf/100)*sgn1;


%Vecteur initial g�n�r� al�atoirement.
Xin=[0, u0, 0, w0, theta0, omega0]; %%[X0;u0;Z0;w0;tetha0;omega0]
%vecteur du temps de simulation
xspan=[0.0:deltat:tmax];


%ODE numerical solution:
  options = odeset('OutputFcn',@odeplot); 
[t, x]=ode45('motion_plate_Shimoi', xspan, Xin),options);


%resultats sous forme de tableaux
Z=zeros(taille,1);
X=zeros(taille,1);

T=zeros(taille,1);

U=zeros(taille,1);
W=zeros(taille,1);

OMega=zeros(taille,1);
THETA=zeros(taille,1);
for i = 1:(taille)
    X(i,1)=x(i,1);
    Z(i,1)=x(i,3);
    U(i,1)=x(i,2);
    W(i,1)=x(i,4);
    OMega(i,1)=x(i,6);
    THETA(i,1)=x(i,5);
    T=t;
end

%%trajectoire moyenne
Xmoy=(Xmoy+X);
Zmoy=(Zmoy+Z);
Umoy=(Umoy+U);
Wmoy=(Wmoy+W);
THetamoy=(THetamoy+THETA);
OMegamoy=(OMegamoy+OMega);



% trajectoire min et max
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

  
%exporter les tableaux dans excel
filename = 'parametriqueVf.xlsx';

% Pour exporter Toutes les it�rations
% indice1=(['A',num2str(j+((j-1)*tmax/deltat))]);
% indice2=(['B',num2str(j+(j-1)*tmax/deltat)]);
% indice3=(['C',num2str(j+(j-1)*tmax/deltat)]);
% indice4=(['D',num2str(j+(j-1)*tmax/deltat)]);
% xlswrite(filename,X,1,['A',num2str(j+((j-1)*tmax/deltat))])
% xlswrite(filename,Z,1,['B',num2str(j+(j-1)*tmax/deltat)])
% xlswrite(filename,T,1,['C',num2str(j+(j-1)*tmax/deltat)])
% xlswrite(filename,theta0,1,['D',num2str(j+(j-1)*tmax/deltat)])


end


xlswrite(filename,Xmoy/n,1,'AH2')
xlswrite(filename,Zmoy/n,1,'AI2')

% xlswrite(filename,Umoy/n,3,'D2')
% xlswrite(filename,Wmoy/n,3,'E2')
% % 
% xlswrite(filename,OMegamoy/n,2,'AK2')
% xlswrite(filename,THetamoy/n,2,'P2')

% xlswrite(filename,Xmax,3,'D2')
% xlswrite(filename,Zmax,3,'E2')
% 
% xlswrite(filename,Xmin,3,'G2')
% xlswrite(filename,Zmin,3,'H2')

% xlswrite(filename,Repartition,3,'R2')

%  xlswrite(filename,Plage,1,'Q2')
toc


end