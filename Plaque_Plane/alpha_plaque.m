function [alpha,beta] = alpha_plaque(u,v,w, vf,vyf,wf)
% Compute the relative velocity of the particule

% first, compute relative velocity component
urel=vf-u;
vyrel=vyf-v;
wrel=wf-w;

vrel=sqrt(urel^2+vyrel^2+wrel^2);
% angles du vent relatif en degrés


alpha =asind(wrel/vrel);
beta = acosd(urel/(vrel*cosd(alpha)+1e-12));


end
