function  [Uq] = velocityTest(Xq, Yq, Zq, component, method) %fonction qui renvoie les composantes de vitesse aux points (Xq, Yq, Zq) connaissant les vitesses 
%aux points du maillage

V1=zeros(100,1)+70;
V2=zeros(100,1);
V3=zeros(100,1);

A=[V1,V2,V3,rx,ry,rz];

P = [A(:,4) A(:, 5) A(:, 6)]  ; % tableau des coordonnées des points du maillage

Mat_Vx =A(:,1); %vecteur colonne contenant les vitesses U en chaque point du maillage

Points_U = scatteredInterpolant(P, Mat_Vx, method); %interpolant de la vitesse U

Mat_Vy =A(:,2); %vecteur colonne contenant les vitesses V en chaque point du maillage

Points_V = scatteredInterpolant(P, Mat_Vy, method); %interpolant de la vitesse V

Mat_Vz =A(:,3); %vecteur colonne contenant les vitesses W en chaque point du maillage

Points_W = scatteredInterpolant(P, Mat_Vz, method); %interpolant de la vitesse W

Vxq = Points_U(Xq, Yq, Zq); %Vitesse U interpolee au point (Xq, Yq, Zq)
Vyq = Points_V(Xq, Yq, Zq); %Vitesse V interpolee au point (Xq, Yq, Zq)
Vzq = Points_W(Xq, Yq, Zq); %Vitesse W interpolee au point (Xq, Yq, Zq)

V= [Vxq;Vyq;Vzq];

if component == 1
   Uq =  V(1)
end
if component == 2 
   Uq = V(2)
end
if component == 3
   Uq = V(3)
end 
end