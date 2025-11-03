function [figure,aw,cmt, Cmf, cm0wing] = UAS1_LongStab(sweep, quarterSweep, wingspan, wingarea, fuselagew,MAC,...
    h, hn, tailsweephalf, tailspan, tailarea,lt)
% UAS 1     Aft-Swept + Tail
% Airfoil Wing Relations
%  Using Values of wing and DAE31 Airfoil

V = 13.889    ; % m/s max speed wanted
a = 343 ; % m/s at ~ 120m (max ceiling)
M = V/a;
Lambdawing = sweep; % sweep
LambdaQwing = quarterSweep; % half sweep angle at mid point line of wing planform
a0wing = 2.006*pi   ; % radians, lift curve slope of airfoil
cd0wing = 0.01331;
cm0wing = -0.1568;
Swing = wingarea    ; % wing area
bwing = wingspan; % wing span
ARwing = bwing^2/Swing;
df = fuselagew      ; % fuselage diamater/width
sw = 1 - 2*(df/bwing)^2     ; % s factor accounts for induced drag caused by span loading
                          % due to fuselage presence
u = 0.99  ; % theoretical oswald factor
Qw = 1/(u*sw);
ew = 1/Qw   ; % Oswald Efficiency factor assuming inviscid flow
Cw = MAC ; % MAC of wing
hw =   0.2 ; % xcg/MAC
hnw = hn   ; % xac of wing/MAC
zw =  0   ; % Z Distance from wing's neutral point and a/c cg 
alpha = linspace(-10, 10, 100000).* pi/180;


% Theoretical Lift Curve of Wing Dependent on Aifoil Lift Curve Slope

kw = 2*pi*(a0wing).^-1;

aw = (pi*ARwing) * (1 + sqrt(1+ ((1-M^2)*(cos(LambdaQwing)))*(((pi*ARwing)/(a0wing*cos(LambdaQwing)))^2))).^-1;
%aw = (2*pi*ARwing) * (2 + sqrt( (ARwing^2*(1-M^2)*(kw)^2)* (1 + (tan(LambdaQwing)^2)/(1-M^2)) + 4)).^-1;
% ^should be lambda at half chord line on wing

clw = aw * alpha;


% Theoretical Drag coefficient using C_d min from a given airfoil

cdw = cd0wing + clw.^2.*(pi * ARwing * ew)^(-1); % full lambda used


% Theoretical Moment Coefficient of wing from cd and cl values

cmw = cm0wing + (clw.*cos(alpha) + cdw.*sin(alpha)).*(hw-hnw) + (clw.*sin(alpha)-cdw.*cos(alpha)).*(zw/Cw);
% 
% figure = plot(alpha.*180/pi,cmw); title("C_m vs \alpha", "UAS 1"); 
% xlabel("angle of attack (^o)"); ylabel("Wing Moment Coefficient");
% grid on; hold off


% Airfoil Tail Relations
% NACA 0012

Lambdahtail = tailsweephalf;
a0tail = 2*pi; % lift-curve slope of airfoil
cd0tail = 0.007;
cm0tail = 0;
Stail = tailarea;
btail = tailspan;
ARtail = btail^2/Stail;
df = fuselagew     ; % fuselage diamater/width
st = 1 - 2*(df/btail)^2     ; % s factor accounts for induced drag caused by span loading
              % due to fuselage presences
u = 0.99  ; % theoretical oswald factor
Qt = 1/(u*st);
et = 1/Qt   ; % Oswald Efficiency factor assuming inviscid flow

% Theoretical Lift Curve of Tail Dependent on Aifoil Lift Curve Slope
kt = 2*pi*(a0tail).^-1;

at = (2*pi*ARtail) * (2 + sqrt( (ARtail^2*(1-M^2)*(kt).^2).* (1 + (tan(Lambdahtail)^2)/(1-M^2)) + 4)).^(-1);

clt = at* alpha;


% Theoretical Drag coefficient using C_d min from a given airfoil

cdt = cd0tail + clt.^2.*(pi * ARtail * et);

% Horizontal Tail Volume Ratio

ltail = lt   ; % Distance from tail aerodynamic center to wing aerodynamic center
Vh_hat = (ltail * Stail)/(Cw*Swing);
Vhtail = Vh_hat - Stail/Swing * (hw - hnw); % Tail effective velocity

% Tail Efficiency

etahtail = (Vhtail/ V);

% Theoretical Moment Coefficient of wing from cd and cl values

cmt = -etahtail * Vh_hat * clt + etahtail * (Stail/Swing) * (hw-hnw) * clt;

% Moment of Whole body

% Moment of Fuselage
EV = 0.020106; % Effective Volume of Fuselage shape
Cmf = 2*EV/(Swing*Cw)*alpha;

% Wing body pitching moment. 
cm0wf = cm0wing * ARwing*cos(Lambdawing)^2/(ARwing + 2*cos(Lambdawing)^2); 
% cm0wing is wing airfoil pitching moment, and full lambda used here

Cm1 = cm0wf + cmw + cmt; % + Cmf; % Full moment coefficient to find longitudinal stability
% plot showing C_m vs alpha; for stability, Cm_alpha is negative
% 
figure = plot(alpha.*180/pi,Cm1); title("C_m vs \alpha", "UAS 1"); 
xlabel("angle of attack (^o)"); ylabel("Total Moment Coefficient");
grid on; hold off
end
