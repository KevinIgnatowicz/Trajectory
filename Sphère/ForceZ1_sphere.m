function  [Fz_sp] = ForceZ1_sphere(Vf,vyf, Wf,rhoair,muair,u,v,w,dsp,alpha)

% Sphere static drag equations
% Calcul de la force selon Z
% g=9.81;
% the drag is parallel to the relative velocity
% the lift is perpendical

% compute relative velocity, 
vrel=sqrt((Vf-u).^2 + (vyf-v).^2 + (Wf-w).^2 );
% compute the relative Reynolds number
Re=rhoair*dsp*vrel/muair;

% compute de drag coefficient and lift
Cd = CD_sp(Re);
Cl = 0.0d0;


% compute the force in absolute z direction
 
% alpha is the wind relative angle of attack
Afront= pi*dsp.^2/4.0;
Fz_sp= 0.5*rhoair*Afront*vrel.^2.0*(Cd*sind(alpha));

end
