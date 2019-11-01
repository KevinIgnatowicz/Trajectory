
variables_plaque;
Srf=xlsread('AileBWB.xlsx');
Pt=[Srf(:,1),Srf(:,2),Srf(:,3)];
sz=size(Pt,1);
global Points_U Points_V Points_W Xin Pr

A = csvread(filename,1,0);

Pts = [A(:,24) A(:, 25) A(:, 26)]  ; % tableau des coordonnées des points du maillage

Mat_Vx =(A(:,2))./(A(:,1)); %vecteur colonne contenant les vitesses U en chaque point du maillage

Points_U = scatteredInterpolant(Pts, Mat_Vx, method); %interpolant de la vitesse U

Mat_Vy =(A(:,3))./(A(:,1)); %vecteur colonne contenant les vitesses V en chaque point du maillage

Points_V = scatteredInterpolant(Pts, Mat_Vy, method); %interpolant de la vitesse V

Mat_Vz =(A(:,4))./(A(:,1)); %vecteur colonne contenant les vitesses W en chaque point du maillage

Points_W = scatteredInterpolant(Pts, Mat_Vz, method); %interpolant de la vitesse W
CI_MAT = zeros (500, 3);
VI_MAT = zeros (500,3);
ANGLESI_MAT = zeros (500,3);
PR_MAT = zeros (500,2);
for i=1:500
    Traj_NB = i
    CI=Pt(randi(sz),:);
    VI=[-0.7*rand,0,sign(CI(3))*0.7*rand];
    anglesI=[(rand-0.5)*randi(180),(rand-0.5)*randi(180),(rand-0.5)*randi(180)];
    [q1I,q2I,q3I,q4I]=Euler2quat(anglesI(1)*pi/180,-anglesI(2)*pi/180,-anglesI(3)*pi/180);
   CI_MAT(i,:) = CI;
   VI_MAT(i,:) = VI;
   ANGLESI_MAT(i,:) = anglesI;
Xin=[0,0,0,q1I,q2I,q3I,q4I,CI(1)+0.001*sign(CI(1)-17),VI(1),CI(2),VI(2),CI(3)+0.001*sign(CI(3)),VI(3),anglesI(1)*pi/180,-anglesI(2)*pi/180,-anglesI(3)*pi/180];
main;
PR_MAT(i,:) = Pr;
end
xlswrite('Classeur2.xlsx', CI_MAT, "Position")
xlswrite('Classeur2.xlsx', VI_MAT, "Vitesse")
xlswrite('Classeur2.xlsx', ANGLESI_MAT, "Angles")
xlswrite('Classeur2.xlsx', PR_MAT, "empreintes")