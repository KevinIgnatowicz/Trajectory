B = csvread('surface_flow.csv',1,0);
figure('Name', 'Trajectoire')
POINTS=[B(:,2), B(:,3), B(:,4)];
pcshow(POINTS,[0 0 0],'MarkerSize',2)
axis equal
title('Trajectoire 3D');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

BWB=alphaShape(B(:,2), B(:,3), B(:,4),2);
figure(1);
hold on
plot(BWB,'Facecolor', 'white');
hold off
%{
figure('Name', 'Print')
k=boundary(B(:,3), B(:,4));
yb=B(:,3);
zb=B(:,4);
plot(yb(k), zb(k))
axis equal
title('Empreintes des trajectoires dans le plan X=1.4m');
xlabel('Y (m)')
ylabel('Z (m)')

figure('Name', 'Print')
k=boundary(B(:,3), B(:,4));
yb=B(:,3);
zb=B(:,4);
plot(yb(k), zb(k))
axis equal
title('Empreintes des trajectoires dans le plan X=2.8m');
xlabel('Y (m)')
ylabel('Z (m)')
%}