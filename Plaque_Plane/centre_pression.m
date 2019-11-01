function [centre_p] = centre_pression(q1,q2,q3,q4,u,v,w,vf,vyf,wf,L,l)%centre_pression(alpha,theta, beta,psi, L, l)

%Cette fonction calcul le centre de pression pour une plaque en 3D selon
%l'article de Richards (Numerical calculation of the three-dimensional
%motion of wind-borne debris, 2008)


%dimensions = [L,l,e] : dimensions de la plaque
%Calcul du centre de pression (x_c ; y_c)
%angles est le vecteur définissant la vitesse relative angles=(alpha, beta)
%alpha = angle d'attaque ; beta = angle d'inclinaison
%*180.0/pi


Ub=Rq(q1,q2,q3,q4)*[vf-u;vyf-v;wf-w];
alpha_local=asind(Ub(3)/sqrt(Ub(1).^2+Ub(2).^2+Ub(3).^2));
beta_local=asind(Ub(2)/(sqrt(Ub(1).^2+Ub(2).^2+Ub(3).^2)*cosd(alpha_local)+1e-12));
%angle relatif zeta permet de calculer ce centre
zeta = atand(l*tand(beta_local)/L);

c=(L*l)/(L*abs(cosd(beta_local))+l*abs(sind(beta_local)));

x_c = -sign(vf-u)*c/4*(90-abs(alpha_local))/90*cosd(zeta);%-sign(vf-u)*
y_c = -sign(vyf-v)*c/4*(90-abs(alpha_local))/90*sind(zeta);%-sign(vyf-v)*

centre_p = [x_c;y_c;0];
end