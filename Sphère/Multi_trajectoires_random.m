variables;
Srf=xlsread('AileBWB.xlsx');
Pt=[Srf(:,1),Srf(:,2),Srf(:,3)];
sz=size(Pt,1);
global Points_U Points_V Points_W Xin

A = textread(filename);

Pts = [A(:,18) A(:, 19) A(:, 20)]  ; % tableau des coordonnées des points du maillage

Mat_Vx =(A(:,2))./(A(:,1)); %vecteur colonne contenant les vitesses U en chaque point du maillage

Points_U = scatteredInterpolant(Pts, Mat_Vx, method); %interpolant de la vitesse U

Mat_Vy =(A(:,3))./(A(:,1)); %vecteur colonne contenant les vitesses V en chaque point du maillage

Points_V = scatteredInterpolant(Pts, Mat_Vy, method); %interpolant de la vitesse V

Mat_Vz =(A(:,4))./(A(:,1)); %vecteur colonne contenant les vitesses W en chaque point du maillage

Points_W = scatteredInterpolant(Pts, Mat_Vz, method); %interpolant de la vitesse W

for i=1:200
    CI=[0,rand,(rand-0.5)*0.04];%Pt(randi(sz),:);
    VI=[velocity(CI(1), CI(2), CI(3), 1),velocity(CI(1), CI(2), CI(3), 2),velocity(CI(1), CI(2), CI(3), 3)];%[-0.5*rand,0,sign(CI(3))*0.7];
   
Xin=[CI(1),VI(1),CI(2),VI(2),CI(3),VI(3),0,0]
TraceXZ;
end