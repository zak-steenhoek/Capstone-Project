function [figure] = Joined_LongStab(wingarea1, wingarea2, wingspan1,wingspan2, MAC1, MAC2,...
    hn1, hn2, aw, aw2, Cmf, Cmfa2, cmt, cmc, cm0wing, cm0wing2)
%% UAS 1 + UAS 2
V = 13.889    ; % m/s max speed wanted
a = 343 ; % m/s at ~ 120m (max ceiling)
M = V/a;
alpha = linspace(-10, 10, 100).* pi/180;
rho = 1.225; % Density at sea level (kg/m^3)
Swing = wingarea1;
Swing2 = wingarea2;
S_n = wingspan1 + wingspan2; % Total Wing Span

C_n = MAC1 + MAC2; % Total MAC
Cw = MAC1; Cw2 = MAC2;

l1 = hn1; % Distance from MAC to cg of Front Wing
l2 = hn2; % Distance from MAC to cg of aft wing
q = 1/2*rho*V^2;
%%%%%%%%% Note: this model does not take into account vertical height
%%%%%%%%% between wing planforms or effects of canards and tails. We
%%%%%%%%% should, however, keep the wings as far apart as possible for
%%%%%%%%% aerodynamic efficiency.

Cm12 = -l1*Swing/(C_n*S_n)*aw*alpha - l2*Swing2*q*aw2/(C_n*S_n)*alpha +...
    Swing*Cw*cm0wing/(S_n*C_n) + q*Swing2*Cw2/(S_n*C_n)*cm0wing2 +...
    Cmf+Cmfa2+cmt+cmc;


% plot showing C_m vs alpha; for stability, Cm_alpha is negative
figure = plot(alpha.*180/pi,Cm12); title("C_m vs \alpha", "Joined UAS"); 
xlabel("angle of attack (^o)"); ylabel("Total Moment Coefficient");
grid on;

end