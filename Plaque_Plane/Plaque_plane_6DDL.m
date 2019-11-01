function [ dx ] = Plaque_plane_6DDL ( t,x )

variables_plaque;
%rotquaternion test ode solution for rotating quaternion
%   Detailed explanation goes here
% x(1,2,3) contains P, Q, R
% rotation velocity in radian P, Q, R 
% x(4,5,6,7) contains q1, q2, q3, q4
% x(1) => P= d(phi)/dt
% x(2) => Q= d(theta)/dt 
% x(3) => R= d(psi)/dt
% x(4) => q1
% x(5) => q2
% x(6) => q3
% x(7) => q4

% x(8) => x
% x(9) => u
% x(10) => y
% x(11) => v
% x(12) => z
% x(13) => w
% x(14) => Phi
% x(15) => theta
 % x(16) => psi
% dx(1,2,3) acceleration angulaire  en x, y, z
% dx(4,5,6,7) derivative of q1, q2, q3, q4

Vf= velocity(x(8),x(10),x(12),1);
Vyf= velocity(x(8),x(10),x(12),2);
Wf = velocity(x(8),x(10),x(12),3);


%[alpha,beta]=alpha_plaque(x(9),x(11),x(13),Vf,Vyf,Wf);


%compute the forces 
F=Forces_plaque_plane(x(4),x(5),x(6),x(7),L,l,rhoair,x(9),x(11),x(13),Vf,Vyf,Wf,Vq);

%compute the moment of the plate 
if choix_moment == 1 
M=Moment_plaque_plane(L,l,x(1),x(2),x(3),x(4),x(5),x(6),x(7),rhoair,x(9),x(11),x(13),Vf,Vyf,Wf,Vq);
end
if choix_moment == 2
M=Moment_plaque_plane_Tachikawa(L,l,e,x(1),x(2),x(3),x(4),x(5),x(6),x(7),rhoair,x(9),x(11),x(13),Vf,Vyf,Wf,Vq);
end

%compute the angular velocity from the acceleration
  
PQR=invIplate*(M-cross(x(1:3,1),Iplate*x(1:3,1)));

%compute integration drift gain correction lambda
lambda=1.-(x(4)^2+x(5)^2+x(6)^2+x(7)^2);

%compute the derivative


dx(1,1)=PQR(1);
dx(2,1)=PQR(2);
dx(3,1)=PQR(3);
dx(4,1)=-0.5*(x(5)*x(1)+x(6)*x(2)+x(7)*x(3))+ lambda*x(4);
dx(5,1)=0.5*(x(4)*x(1)+x(6)*x(3)-x(7)*x(2))+ lambda*x(5);
dx(6,1)=0.5*(x(4)*x(2)+x(7)*x(1)-x(5)*x(3))+ lambda*x(6);
dx(7,1)=0.5*(x(4)*x(3)+x(5)*x(2)-x(6)*x(1))+ lambda*x(7);

%relations quaternions/angles d'Euler/vitesses de rotation
dx(14,1)=x(1);
dx(15,1)=x(2);
dx(16,1)=x(3);
%[x(16), x(15), x(14)]= quat2euler(x(4), x(5), x(6), x(7));
[x(4), x(5), x(6), x(7)]=Euler2quat(x(14), -x(15), -x(16));

%trajectoires
dx(8,1)=x(9);
dx(9,1)=F(1)/mass;

dx(10,1)=x(11);
dx(11,1)=F(2)/mass;

dx(12,1)=x(13);
dx(13,1)=F(3)/mass-g;

end

