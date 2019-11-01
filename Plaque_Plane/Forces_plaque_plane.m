function [Fi] = Forces_plaque_plane(q1,q2,q3,q4,L,l,rhoair,u,v,w,vf, vyf, wf,Vq)

%vitesse relative dans le référentiel mobile de la plaque

Ub=Rq(q1,q2,q3,q4)*[vf-u;vyf-v;wf-w];

%forces dans le référentiel mobile de la plaque
%hypothèse  : seule la force normale est considérée comme prépondérante sur
%les autres


Fb = 0.5*rhoair*L*l*(Ub(1).^2+Ub(2).^2+Ub(3).^2).*[0;0;CN_3D(q1,q2,q3,q4,u,v,w,vf,vyf,wf,Vq)]; 

%Forces dans le référentiel fixe :

Fi=transpose(Rq(q1,q2,q3,q4))*Fb;

end
