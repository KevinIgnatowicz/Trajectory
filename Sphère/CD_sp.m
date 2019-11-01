function out = CD_sp(Re)
% Calcul du coefficient de train�e de la sph�re
% IN : Re
% OUT : CD (different from previous CD_sp function
% D'apres Bubbles, Drops, and Particle p112 table 5.2 a,b par Clift Grace &
% Weber (1978)



if ( Re < 0.01 )
    
    Cd = 3.0/16.0 + 24.0./(Re+1e-12);
    out = Cd;
   return
    
end

if ( Re<=20 )
    
    w = log10(Re+1e-12);
    pow=(0.82-0.05.*w);
    Cd = (24.*(1+0.1315*Re.^pow))/(Re+1e-12);
    out = Cd;
    return
    
end

if  ( Re <= 260 )
    
    Cd = (24.*(1+0.1935*Re.^(0.6305)))/(Re+1e-12);
   out = Cd;
    return
    
end

if  ( Re <= 1500 )
    
    w = log10(Re+1e-12);
    Cd = 10.^(1.6435-1.1242.*w+0.1558.*w.^2);
    out = Cd;
    return
    
end

if ( Re <= 1.2e4 )
    
    w = log10(Re+1e-12);
    Cd = 10.^(-2.4571+2.5558.*w-0.9295.*w.^2+0.1049.*w.^3);
    out = Cd;
    return
    
end

if ( Re <= 4.4e4 )
    
    w = log10(Re+1e-12);
    Cd = 10.^(-1.9181+0.6370.*w-0.0636.*w.^2);
    out = Cd;
    return
    
end

if( Re <= 3.38e5 )
    
    w = log10(Re+1e-12);
    Cd = 10.^(-4.3390+1.5809.*w-0.1546.*w.^2);
    out = Cd;
    return
    
end

if  ( Re <= 4e5 )
    
    w = log10(Re+1e-12);
    Cd = 29.78 - 5.3.*w;
    out = Cd;
    return
    
end

if( Re <= 1e6 )
    
    w = log10(Re+1e-12);
    Cd = 0.1.*w-0.49;
    out = Cd;
    return
    
end

if ( Re > 1e6 )
    
    Cd = 0.19-8e4./(Re+1e-12);
    out = Cd;
    return
    
end


end

