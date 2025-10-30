%% UAS 2 Forward Swept + Canard
function [figure] = UAS2_LongStab(sweep, quarterSweep, wingspan, wingarea, fuselagew,MAC,...
    h, hn, tailsweephalf, tailspan, tailarea,lt)

% Airfoil Wing Relations
Lambdawing2 = sweep;
LambdaQwing2 = quarterSweep;
a0wing2 = 6    ; % radians
cd0wing2 = 0.00533;
cm0wing2 = -0.0787;
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
hw2 = h ; % MAC/xcg
hnw2 = hn    ; % MAC/xac of wing

% Theoretical Lift Curve of Wing Dependent on Aifoil Lift Curve Slope

kw2 = 2*pi*(a0wing2).^-1;

%aw2 = (2*pi*ARwing2) * (2 + sqrt( (ARwing2^2*(1-M^2)*(kw2).^2).* (1 + (tan(Lambdawing2)^2)/(1-M^2)) + 4)).^(-1);
aw2 = (pi*ARwing2) * (1 + sqrt( ((1-M^2)*(a0wing2*cos(LambdaQwing2)^2))* (1 - ((pi*ARwing2)/(cos(LambdaQwing2)^2))^2))).^-1;

clw2 = aw2* alpha;


% Theoretical Drag coefficient using C_d min from a given airfoil

cdw2 = cd0wing2 + clw2.^2./(pi * ARwing2 * ew2);


% Theoretical Moment Coefficient of wing from cd and cl values

cmw2 = -cm0wing2 - (clw2.*cos(alpha) + cdw2.*sin(alpha)).*(hw2+hnw2) - (clw2.*sin(alpha)-cdw2.*cos(alpha)).*(zw2/Cw2);


% Airfoil Tail Relations %

Lambdahcanard = tailsweephalf;
a0canard = 2*pi;
cl0canard = 0; 
cd0canard = 0.007;
cm0canard = 0;
Scanard = tailspan;
bcanard = tailarea;
ARcanard = bcanard^2/Scanard;
df = fuselagew   ; % fuselage diamater/width
sc = 1 - 2*(df/bcanard)^2     ; % s factor accounts for induced drag caused by span loading
                                % due to fuselage presence
u = 0.99  ; % theoretical oswald factor
Qc = 1/(u*sc);
ec = 1/Qc   ; % Oswald Efficiency factor assuming inviscid flow

% Theoretical Lift Curve of Tail Dependent on Aifoil Lift Curve Slope
kc = 2*pi*(a0canard).^-1;

ac = (2*pi*ARcanard) * (2 + sqrt( (ARcanard^2*(1-M^2)*(kt).^2).* (1 + (tan(Lambdahcanard)^2)/(1-M^2)) + 4)).^(-1);

clc = ac* alpha;


% Theoretical Drag coefficient using C_d min from a given airfoil

cdc = cd0canard + clc.^2./(pi * ARcanard * ec);

% Canard Volume Ratio

lcanard = lt; % Distance from canard aerodynamic center to A/C cg
Vh_hatc = (lcanard * Scanard)/(Cw2*Swing2);


% Canard moment %

cmc = Vh_hatc*cl0canard + Vh_hatc*ac;



% Moment of whole body %

% Moment of Fuselage
EV2 = 6.28319*10^-3; % Effective Volume of Fuselage shape
Cmfa2 = 2*EV2/(Swing*Cw); % Fuselage body moment angle of attack slope (destabilizing)

% Wing body pitching moment. 
cm0wf2 = cm0wing2 * ARwing2*cos(Lambdawing2)^2/(ARwing2 + 2*cos(Lambdawing2)^2); % c_m0 is wing airfoil pitching moment

Cm2 = cm0wf2 + cmw2 + cmc +  Cmfa2;


% plot showing C_m vs alpha; for stability, Cm_alpha is negative
figure = plot(alpha.*180/pi,Cm2); title("C_m vs \alpha", "UAS 2"); 
xlabel("angle of attack (^o)"); ylabel("Total Moment Coefficient");
grid on; hold off;
end