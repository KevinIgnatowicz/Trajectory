angles_ini=[0 5 0; 10 12 0; 20 3 20; 0 10 0; 15 5 10;0 4 0; 14 2 5;10 10 0;15 5 20];
%CI =[P0 Q0 R0 X0 u0 Y0 v0 Z0 w0]
CI= [0 0 0 16.57 -0.5 11.73 0 0.21 0.5;0 0 0 15.87 -0.5 10.96 0 0.231 0.5;0 0 0 15.03 -0.5 10.85 0 -0.22 -0.5;0 0 0 11.87 -0.5 7.61 0 0.26 0.5;0 0 0 18.21 -0.5 13.89 0 0.18 0.5;0 0 0 12.93 -0.5 4.67 0 1.233 0.5;0 0 0 16.16 0 10.01 0 0.22 0.5;0 0 0 15.21 -0.5 11.63 0 -0.14 -0.5;0 0 0 7.45 -0.5 4.38 0 -0.68 -0.5]

%global Xin
for i=1:size(angles_ini,1)
   % [q10, q20, q30, q40]=Euler2quat(angles_ini(i,1),-angles_ini(i,2),-angles_ini(i,3));
   % Xin=[CI(i,1),CI(i,2),CI(i,3),q10,q20,q30,q40,CI(i,4),CI(i,5),CI(i,6),CI(i,7),CI(i,8),CI(i,9),angles_ini(i,1)*pi/180,-angles_ini(i,2)*pi/180,-angles_ini(i,3)*pi/180]
    main;
end
