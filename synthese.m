close all
alpha1 = 1;
alpha2 = 100;
sig0 = 1;
sig1 = 100;
Te = 0.3;
s = tf('s');
%% modélisation du système
A = [0 0 1; 0 -alpha1 0; 0 sig1 -alpha2];
B = [0; sig0; 0];
E = [0; 0; alpha2];
C = [1 0 0];
D = 0;
sysbo = ss(A, [B E], C, [D zeros(size(C,1), size(E,2))]);
sysbod = c2d(sysbo, Te);
Gs = 1/(s*(s+1));
Gz = c2d(Gs,Te);
Ad = sysbod.A; Bd = sysbod.B; Cd = sysbod.C; Dd = sysbod.D;

%% RST

trdes = 5;
xsides = 0.4; % D = 25%
wndes = 2*pi*1.4/trdes;
taudes = trdes/3;
waux = 10*wndes;
Pdes1 = [1/wndes^2 2*xsides/wndes 1];
Pdes2 = conv([1/waux 1], [1/waux 1]);
Pdes = conv(Pdes1, Pdes2);
polescontinus = roots(Pdes);
polesdiscrets = exp(polescontinus*Te);
Pdez = poly(polesdiscrets);
global p1 p2 p3 p4
p1 = Pdez(2);
p2 = Pdez(3);
p3 = Pdez(4);
p4 = Pdez(5);
aide_resolution_eq_bezout;
sim('sythese_sim.slx')
figure()
hold on
plot(hout.time, hout.signals.values); title('h en fonction du temps RST');
yline(1.05, 'r--');
yline(0.95, 'r--');
yline(1.25, 'black');
xline(5, 'black--');
hold off


%% commande par retour d'état

Pol_car = (s + 100)*(s^2 + 2*xsides*wndes*s + wndes^2);
P = zero(Pol_car);
F = acker(A, B, P);
invABF = inv(- A + B*F);
G = 1/(C * invABF * B);
sysbf = ss(A-B*F, [B*G], C, 0);
Pol_car2 = (s + 101)*(s + 100)*(s^2 + 2*xsides*wndes*s + wndes^2);
P2 = zero(Pol_car2);
L = acker(A', C',P)';
Ae = [A zeros(3,1);-C 0];
Be = [B; 0];
Fe = place(Ae, Be, P2);
Fobs = Fe(1:3);
Hobs = -Fe(4);



figure();
hold on
% xlim([0 Tbf(end)])
step(sysbf); title('réponse à un échelon unitaire - commande par retour d état');
yline(1.05, 'r--');
yline(0.95, 'r--');
yline(1.25, 'black');
xline(5, 'black--');
hold off

Pd = exp(P*Te);
Fd = acker(Ad, Bd(:,1), Pd);
invABFd = inv(eye(3) - Ad + Bd(:,1)*Fd);
Gd = 1/(Cd * invABFd * Bd(:,1));
sysbfd = ss(Ad-Bd(:,1)*Fd, [Bd(:,1)*Gd], Cd, 0, Te);
P2d = exp(P2*Te);
Ld = acker(Ad', Cd', Pd)'; %pour discret
Aed = [Ad zeros(3,1);-Cd 1];
Bed = [Bd(:,1); 0];
Fed = acker(Aed, Bed, P2d);
Fobsd = Fed(1:3);
Hobsd = -Fed(4);


figure();
hold on
% xlim([0 Tbfd(end)])
step(sysbfd); title('réponse à un échelon unitaire - commande par retour d état - en discret');
% stairs(Tbfd, hbfd); title('réponse à un échelon unitaire - commande par retour d état - en discret');
yline(1.05, 'r--');
yline(0.95, 'r--');
yline(1.25, 'black');
xline(5, 'black--');
hold off




sysobs = ss(A-L*C-B*F, [B*G], [C], [0]);
sim('synthese_sim_retour_etat.slx', 10)

figure()
hold on
plot(hre.time, hre.signals.values); title('h en fonction du temps retour état');
yline(1.05, 'r--');
yline(0.95, 'r--');
yline(1.25, 'black');
xline(5, 'black--');
hold off