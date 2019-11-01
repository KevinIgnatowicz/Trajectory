function [Mb]=Moment_plaque_plane( L,l,P,Q,R,q1,q2,q3,q4,rhoair,u,v,w,vf,vyf,wf,Vq)


%vitesse dans le repère de la plaque
Ub=Rq(q1,q2,q3,q4)*[vf-u;vyf-v;wf-w];

%centre de pression :
C_p=centre_pression(q1,q2,q3,q4,u,v,w,vf,vyf,wf,L,l);

%Force apliquée au centre de pression : 
Fb=Rq(q1,q2,q3,q4)*Forces_plaque_plane(q1,q2,q3,q4,L,l,rhoair,u,v,w,vf,vyf,wf,Vq);

%Moment de la force appliquée au centre de pression :
Mf=cross(C_p,Fb);

%Damping moment (moment d'amortissement)
Cdm=[-0.2;-0.185;-0.01];

Md=Cdm.*0.5*rhoair*(sqrt(Ub(1).^2+Ub(2).^2+Ub(3).^2)+sqrt(P.^2+Q.^2+R.^2)*l/2)*L*l.^3.*[P;Q;R];

%Moment total

Mb=Mf+Md;

end



