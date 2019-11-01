function [alpha,betha] = alpha(u,v, w, vf,vyf, wf)
% Compute the relative velocity of the particule

% first, compute relative velocity components
urel=vf-u;
vr= vyf-v; 
wrel=wf-w;
%relative velocity magnitude
vrel = sqrt(urel.^2+vr.^2+wrel.^2)

% angles du vent relatif en degres

alpha = asind(wrel/vrel); %angle d'attaque par rapport au plan XY
betha = asind(vr/(vrel*cosd(alpha)+1e-12));%angle d'inclinaison, par rapport au plan XZ


end
