function [spline_Cn]=CN_3D(q1,q2,q3,q4,u,v,w,vf,vyf,wf,Vq)

Ub=Rq(q1,q2,q3,q4)*[vf-u;vyf-v;wf-w];
alpha_local=asind(Ub(3)/(sqrt(Ub(1).^2+Ub(2).^2+Ub(3).^2)+1e-12));
beta_local=asind(Ub(2)/(sqrt(Ub(1).^2+Ub(2).^2+Ub(3).^2)*cosd(alpha_local)+1e-12));

if  alpha_local>=0
spline_Cn = Vq(alpha_local, abs(beta_local));

end
if  alpha_local<0
    spline_Cn = -Vq(-alpha_local,abs(beta_local));
   
end


end