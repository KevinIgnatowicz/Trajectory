function  [CN_theta] = CN_Holmes(gamma)

% Holmes equations
% Calculation of the normal force coefficient in terms of the angle of
% attack 
% gamma angle d'attaque total


 for i=1:length(gamma)
    
if (mod(gamma(i),180)<40)                                %%When theta<=40�
    
    CN_theta(i) = 1.7*(mod(gamma(i),180)./40)    ;       %%Eq.39a      
    
elseif (40 <= mod(gamma(i),180)) && (mod(gamma(i),180) <=140)     %%When 40�<=theta<=140�                        
   
     CN_theta(i) =1.15;                         %%Eq.39b 
   
else                                            %%When 140�<theta<=180�
    
  CN_theta(i) = 1.7*(180-mod(gamma(i),180))./40  ;      %%Eq.39c

end
 end

% a='CN_Holmes'
end
