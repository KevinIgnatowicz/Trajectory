function [ dx ] = rotquaternion( t,x )
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
% dx(1,2,3) acceleration angulaire  en x, y, z
% dx(4,5,6,7) derivative of q1, q2, q3, q4



%compute the angular velocity from the acceleration
global Iplate invIplate  
PQR=invIplate*(M-cross(x(1:3,1),Iplate*x(1:3,1)));


%compute integration drift gain correction lambda
lambda=1.-(x(4)^2+x(5)^2+x(6)^2+x(7)^2);

%compute the derivative 
P=x(1) % not needed but to clarify the derivative formula
Q=x(2) %
R=x(3) %
dx(1,1)=PQR(1);
dx(2,1)=PQR(2);
dx(3,1)=PQR(3);
dx(4,1)=-0.5*(x(5)*P+x(6)*Q+x(7)*R)+ lambda*x(4);
dx(5,1)=0.5*(x(4)*P+x(6)*R-x(7)*Q)+ lambda*x(5);
dx(6,1)=0.5*(x(4)*Q+x(7)*P-x(5)*R)+ lambda*x(6);
dx(7,1)=0.5*(x(4)*R+x(5)*Q-x(6)*P)+ lambda*x(7);

    

end

