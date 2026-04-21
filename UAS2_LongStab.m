%% UAS 2 Forward Swept + Canard
function [fig, aw2, cmc, cmt, Cmf2, cm0wing2, cl_tot2, cd_tot2] = UAS2_LongStab(sweep, quarterSweep, wingspan, wingarea, fuselagew,MAC,...
    hnw, hnwb, tailsweephalf, tailquartersweep, tailspan, tailarea,lt)

V = 13.889    ; % m/s max speed wanted
a = 343 ; % m/s at ~ 120m (max ceiling)
M = V/a;
alpha = linspace(-10, 10, 100000).* pi/180 ; %- alpha;
alpha1 = alpha -6.82 * pi/180;
% Airfoil Wing Relations
Lambdawing2 = sweep;
LambdaQwing2 = quarterSweep;
a0wing2 = 6;% 2.006*pi    ; % radians (6)
cd0wing2 = 0.01331; % (0.00533)
cm0wing2 = -0.1568; % (-0.0787)
cl0 = 0.7288;
Swing2 = wingarea;
bwing2 = wingspan;
ARwing2 = bwing2^2/Swing2;
df2 = fuselagew  ; % fuselage diamater/width
sw2 = 1 - 2*(df2/bwing2)^2     ; % s factor accounts for induced drag caused by span loading
              % due to fuselage presence
u = 0.99  ; % theoretical oswald factor
Qw2 = 1/(u*sw2);
ew2 = 1/Qw2   ; % Oswald Efficiency factor assuming inviscid flow
zw2 = 0  ;   % Z Distance from wing's neutral point and a/c cg 
Cw2 = MAC  ; % MAC of wing
hw2 = 0.20 ; % xcg/MAC
hnw2 = hnw/MAC    ; % xac of wing/MAC

% Theoretical Lift Curve of Wing Dependent on Aifoil Lift Curve Slope

%aw2 = (2*pi*ARwing2) * (2 + sqrt( (ARwing2^2*(1-M^2)*(kw2).^2).* (1 + (tan(Lambdawing2)^2)/(1-M^2)) + 4)).^(-1);
aw2 = (pi*ARwing2) * (1 + sqrt( ((1-M^2)*(a0wing2*cos(LambdaQwing2)^2))* (1 - ((pi*ARwing2)/(cos(LambdaQwing2)^2))^2))).^-1;
aw2 = real(aw2);
clw2 = cl0 + real(aw2* alpha1);


% Theoretical Drag coefficient using C_d min from a given airfoil

cdw2 = cd0wing2 + real(clw2.^2./(pi * ARwing2 * ew2));


% Theoretical Moment Coefficient of wing from cd and cl values

cmw2 = +cm0wing2 + (clw2.*cos(alpha1) + cdw2.*sin(alpha1)).*(hw2+hnw2) + (clw2.*sin(alpha1)-cdw2.*cos(alpha1)).*(zw2/Cw2);

% % plot showing C_m vs alpha; for stability, Cm_alpha is negative
% figure = plot(alpha.*180/pi,cdw2); title("C_d vs \alpha", "UAS 2"); 
% xlabel("angle of attack (^o)"); ylabel("Drag Coefficient");
% grid on; hold off;

% Airfoil Tail Relations %

Lambdahcanard = tailsweephalf;
LambdaQcanard = tailquartersweep;
a0canard = 2*pi;
cl0canard = 0; 
cd0canard = 0.007;
cm0canard = 0;
bcanard = tailspan;
Scanard = tailarea;
ARcanard = bcanard^2/Scanard;
df = fuselagew   ; % fuselage diamater/width
sc = 1 - 2*(df/bcanard)^2     ; % s factor accounts for induced drag caused by span loading
                                % due to fuselage presence
u = 0.99  ; % theoretical oswald factor
Qc = 1/(u*sc);
ec = 1/Qc   ; % Oswald Efficiency factor assuming inviscid flow

% Theoretical Lift Curve of Tail Dependent on Aifoil Lift Curve Slope

% ac = (2*pi*ARcanard) * (2 + sqrt( (ARcanard^2*(1-M^2)*(kc).^2).* (1 + (tan(Lambdahcanard)^2)/(1-M^2)) + 4)).^(-1);
ac = (pi*ARcanard) * (1 + sqrt( ((1-M^2)*(2*pi*cos(LambdaQcanard)^2))* (1 - ((pi*ARcanard)/(cos(LambdaQcanard)^2))^2))).^-1;

clc = ac* alpha;


% Theoretical Drag coefficient using C_d min from a given airfoil

cdc = cd0canard + clc.^2./(pi * ARcanard * ec);

% Canard Volume Ratio

lcanard = lt; % Distance from canard aerodynamic center to A/C cg
Vh_hatc = (lcanard * Scanard)/(Cw2*Swing2);


% Canard moment %

cmc = Vh_hatc*cl0canard + Vh_hatc*ac;


LambdaQtail = tailsweephalf;
a0tail = 2*pi; % lift-curve slope of airfoil
cd0tail = 0.007; %0.007
cm0tail = 0; %0
%cd0tail = 0.01331;
%cm0tail = -0.1568;
Stail = tailarea;
btail = tailspan;
ARtail = btail^2/Stail;
df = fuselagew     ; % fuselage diamater/width
st = 1 - 2*(df/btail)^2     ; % s factor accounts for induced drag caused by span loading
              % due to fuselage presences
u = 0.99  ; % theoretical oswald factor
Qt = 1/(u*st);
et = 1/Qt   ; % Oswald Efficiency factor assuming inviscid flow

% Downwash angle of tail
epsilon = 2 *alpha.* aw2* (pi*ARwing2)^(-1); % d(epsilon)/d(alpha)
epsilon0 = 2*pi*aw2*(- 6.82 * pi/180)*(pi*ARwing2).^(-1); % downwash angle at zero angle of attack
it = 0 * pi/180; % Angle of incidence -- NOTE: 

% Theoretical Lift Curve of Tail Dependent on Aifoil Lift Curve Slope

%at = (2*pi*ARtail) * (2 + sqrt( (ARtail^2*(1-M^2)*(kt).^2).* (1 + (tan(Lambdahtail)^2)/(1-M^2)) + 4)).^(-1);
at = (pi*ARtail) * (1 + sqrt(1+ ...
    ((1-M^2)*(cos(LambdaQtail)))*(((pi*ARtail)/...
    (a0tail*cos(LambdaQtail)))^2))).^-1;

clt = at * alpha.*(1 - epsilon) - at*(epsilon0 + it);



% Theoretical Drag coefficient using C_d min from a given airfoil

cdt = cd0tail + clt.^2./(pi * ARtail * et);

% Horizontal Tail Volume Ratio

ltail =  0.07   ; % Distance from tail aerodynamic center to wing aerodynamic center
Vh_hat = (ltail * Stail)/(Cw2*Swing2);
Vhtail = Vh_hat - Stail/Swing2 * (hw2 - hnw2); % Tail effective volume

% Tail Efficiency

etahtail = 0.9;


% Theoretical Moment Coefficient of wing from cd and cl values

cmt = -etahtail * Vh_hat * clt + etahtail * (Stail/Swing2) * (hw2-hnw2) * clt;
% Moment of whole body %

% Moment of Fuselage
EV2 = 6.28319*10^-3; % Effective Volume of Fuselage shape
Cmf2 = 2*EV2/(Swing2*Cw2)*alpha; % Fuselage body moment angle of attack slope (destabilizing)

% Wing body pitching moment. 
cm0wf2 = cm0wing2 * ARwing2*cos(Lambdawing2)^2/(ARwing2 + 2*cos(Lambdawing2)^2); % c_m0 is wing airfoil pitching moment

Cm2 = cm0wf2 - cmw2 + cmc +cmt + Cmf2;

% lift and drag totals
cl_tot2 = clw2 + clc +clt;
cd_tot2 = cdw2 + cdc +cdt;


% plot showing C_m vs alpha; for stability, Cm_alpha is negative
% figure = plot(alpha.*180/pi,Cm2); title("C_m vs \alpha", "UAS 2"); 
% xlabel("angle of attack (^o)"); ylabel("Total Moment Coefficient");
% grid on; hold off;

fig = figure()
fontsize(gcf, scale = 10)
set(0,'DefaultAxesFontName','Times New Roman')
cmytot3 = [0.13 0.09 0.05 0.01 -0.02 -0.05 -0.09 -0.12 -0.15 -0.18 -0.2 -0.225 -0.25 -0.275 -0.3 -0.31 -0.33 -0.35 -0.36 -0.38];
alpha2 = linspace(-10, 10, 20);
plot(alpha2, cmytot3)
hold on; grid on; set(gcf,'Color','w');
plot(alpha.*180/pi,Cm2); 
xlabel('Angle of Attack (deg)')
ylabel('Total Pitching Moment Coefficient')
title('Payload Aircraft Longitudinal Stability')
% ylim([-1 1]); 
xlim([-10 10])
legend('OpenVSP Data', 'Calculated Data')

hold off;

end