function [Rq] = Rq(q1,q2,q3,q4)

Rq=zeros(3,3);

Rq(1)= q1^2+q2^2-q3^2-q4^2;
Rq(2)= 2*q2*q3-2*q1*q4;
Rq(3)= 2*q2*q4+2*q1*q3;
Rq(4)= 2*q2*q3+2*q1*q4;
Rq(5)= q1^2-q2^2+q3^2-q4^2;
Rq(6)= 2*q3*q4-2*q1*q2;
Rq(7)= 2*q2*q4-2*q1*q3;
Rq(8)= 2*q3*q4+2*q1*q2;
Rq(9)= q1^2-q2^2-q3^2+q4^2;

end