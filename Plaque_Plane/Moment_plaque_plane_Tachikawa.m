function [Mb]=Moment_plaque_plane_Tachikawa( L,l,e,P,Q,R,q1,q2,q3,q4,rhoair,u,v,w,vf,vyf,wf,Vq)


%vitesse dans le repère de la plaque
Ub=Rq(q1,q2,q3,q4)*[vf-u;vyf-v;wf-w];

%centre de pression :
C_p=centre_pression(q1,q2,q3,q4,u,v,w,vf,vyf,wf,L,l);

%Force apliquée au centre de pression : 
Fb=Rq(q1,q2,q3,q4)*Forces_plaque_plane(q1,q2,q3,q4,L,l,rhoair,u,v,w,vf,vyf,wf,Vq);

%Moment de la force appliquée au centre de pression :
Mf=cross(C_p,Fb);

%Damping moment (moment d'amortissement)
if L/l == 1
W0=0.64*sqrt(Ub(1)^2+Ub(2)^2+Ub(3)^2)/L; %W0 = 0.32*2V/L plaque carrée
W1=0.64*sqrt(Ub(1)^2+Ub(2)^2+Ub(3)^2)/L;
end
if L/l == 2
W0=(-1.25*(e/l)+0.5163)*2*sqrt(Ub(1)^2+Ub(2)^2+Ub(3)^2)/l; %W0 = 0.48*2V/l plaque rectangulaire
W1=0.8*sqrt(Ub(1)^2+Ub(2)^2+Ub(3)^2)/L;
end


LSCmr = [l*L*l*sign(P)*Tachikawa(abs(P)/W0,L,l);L*L*l*sign(Q)*Tachikawa(abs(Q)/W1,L,l);L*L*e*sign(R)*Tachikawa(abs(R)/W0,L,l)];
Md = 0.5*rhoair*(Ub(1)^2+Ub(2)^2+Ub(3)^2)*LSCmr;

%Moment total

Mb=Mf+Md;
end
