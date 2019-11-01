function [Vq]=initialisation_Cn(SR)
%On initialise le coefficient Cn afin de pas le refaire à chaque itérations
%pour gagner en temps de calcul

%On commence par importer les données correspondant au ratio SR donné en
%entrée
if SR == 1
    M = importdata('sr1.txt');
elseif SR==2
    M = importdata('sr2.txt');
elseif SR == 4
    M = importdata('sr4.txt');
end

%La première et troisième colonne correspondent aux angles d'attaque et
%d'inclinaison
alp = M(:,1);
bet = M(:,3);


% length(bet);
% n=length(alp);



Cn = M(:,2);

size(Cn);

Vq = scatteredInterpolant(alp,bet, Cn);



end