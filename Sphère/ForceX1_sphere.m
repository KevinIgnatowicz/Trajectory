function  [Fx_sp] = ForceX1_sphere(Vf,vyf, Wf,rhoair,muair,u,v,w,dsp,alpha,betha)

% Sphere static drag equations
% Calcul de la force selon X
% the drag is parallel to the relative velocity
% the lift is perpendical

% compute relative velocity, 
vrel=sqrt((Vf-u).^2 + (vyf-v).^2 + (Wf-w).^2 );
% compute the relative Reynolds number
Re=rhoair*dsp*vrel/muair;

% compute de drag coefficient and lift
Cd = CD_sp(Re);
Cl = 0.0d0;


% compute the force in absolute x direction
% for generality, wind relatie angle is taken into account 
% alpha and betha are the wind relative angles
Afront= pi*dsp.^2/4.0;
if Vf-u==0
    a=1;
else a=sign(Vf-u);
end
Fx_sp= a*0.5*rhoair*Afront*vrel.^2.0*(Cd.*cosd(betha)*cosd(alpha));

end
