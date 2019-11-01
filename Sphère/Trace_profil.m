A = csvread('Profil_NACA0012.csv')

figure('Name', 'Trajectoire')
plot(A(:,1), A(:,2))
axis equal