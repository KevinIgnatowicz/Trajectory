function  [Uq] = velocity(Xq, Yq, Zq, component) %fonction qui renvoie les composantes de vitesse aux points (Xq, Yq, Zq) connaissant les vitesses 
%aux points du maillage

%global Points_U Points_V Points_W

%Vxq = Points_U(Xq, Yq, Zq); %Vitesse U interpolee au point (Xq, Yq, Zq)
%Vyq = Points_V(Xq, Yq, Zq); %Vitesse V interpolee au point (Xq, Yq, Zq)
%Vzq = Points_W(Xq, Yq, Zq); %Vitesse W interpolee au point (Xq, Yq, Zq)

%V= [Vxq;Vyq;Vzq];

if component == 1
   Uq = 0.0000000000001% V(1)
end
if component == 2 
   Uq = 0.0000000000001% V(2)
end
if component == 3
   Uq = -1% V(3)
end 
end
