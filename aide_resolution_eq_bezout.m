Te = 0.3;
trdes = 5;

%%
% Définition des paramètres du modèle 
tau=1;
global p1 p2 p3 p4
a1= -1 - exp(-Te/tau) ;
a2= exp(-Te/tau);
b1= Te - tau + tau*exp(-Te/tau);
b2= tau - exp(-Te/tau)*(tau + Te);

%% 
% Définition des performances désirées 
%%
% Résolution manuelle de l'équation de Bezout
ap1=a1-1;
ap2=a2-a1;
ap3=-a2;

Pdes=[1 p1 p2 p3 p4];
x=inv([1 0 0 0 0; ap1 1 b1 0 0; ap2 ap1 b2 b1 0;ap3 ap2 0 b2 b1;...
   0 ap3 0 0 b2])*Pdes';
% paramètres du correcteur
global sp1 r0 r1 r2
sp1=x(2); % il s'agit de s'1 du cours.
r0=x(3);
r1=x(4);
r2=x(5);