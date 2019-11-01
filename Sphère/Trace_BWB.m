B = csvread('Surface_BWB.csv',1,0);
figure('Name', 'Trajectoire')
POINTS=[B(:,24), B(:,25), B(:,26)];
pcshow(POINTS,[0 0 0],'MarkerSize',10)
axis equal
title('Trajectoire 3D');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

figure('Name', 'Print')
k=boundary(B(:,25), B(:,26));
yb=B(:,25);
zb=B(:,26);
plot(yb(k), zb(k))
axis equal
title('Emprintes des trajectoires dans le plan X=20m');
xlabel('Y (m)')
ylabel('Z (m)')

