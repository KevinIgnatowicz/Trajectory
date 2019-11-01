function [x, y, z]=quat2euler(q1,q2,q3,q4)
% convert quaternion angle to euler angle
        %compute theta in radian
        y= asin(-2.*(q2.*q4-q1.*q3));
        %compute psi in radian
        x= acos((q1.^2+q2.^2-q3.^2-q4.^2)./cos(y)).*sign(2.*(q2.*q3+q1.*q4));
        %compute phi
        z= acos((q1.^2-q2.^2-q3.^2+q4.^2)./cos(y)).*sign(2.*(q3.*q4+q1.*q2));
end

