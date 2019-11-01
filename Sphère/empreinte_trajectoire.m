function [Print]= empreinte_trajectoire(X,Y,Z,X_plan)
%Fonction qui renvoie les coordonnées (Y_Print,Z_Print) du point d'intersection entre
%la trajectoire définie par les points contenus dans [X Y Z] et le plan de
%normale x positionné à l'abscisse X_plan

A=find(X<X_plan);

if size(A,1)==size(X,1)
    Print='NaN'
elseif size(A,1)==0
    Print='NaN'
else
B=max(A);

Y_Print = ((Y(B+1)-Y(B))/(X(B+1)-X(B)))*(X_plan-X(B))+Y(B);
Z_Print = ((Z(B+1)-Z(B))/(X(B+1)-X(B)))*(X_plan-X(B))+Z(B);

Print = [Y_Print, Z_Print]
end
end
