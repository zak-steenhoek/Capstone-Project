%% Aricraft Nuetral Points
clc; clear; close all;

% UAS 1
width = 0.150; df = width;
length = 1.000;
UAS1x = linspace(0,1,100);
t = width/length;
UAS1Fuselage = 5 * t * (0.2969 * sqrt(UAS1x) - 0.1260 * UAS1x - 0.3516 * UAS1x.^2 + 0.2843 * UAS1x.^3 - 0.1015 * UAS1x.^4);
UAS1x = [UAS1x fliplr(UAS1x)] .* length;
UAS1Fuselage = [UAS1Fuselage -fliplr(UAS1Fuselage)];

% Wing UAS 1
NoseSetbackDist = 0.200; % [m]
Root_Chord = 0.2485; % [m] 0.240
Tip_Chord = 0.1175; % [m] 0.060
Half_Span = 1/2; % [m]
wspan1 = 2 * Half_Span;
Sweep_Angle = 26.2489 * (pi/180); Lambdaw1 = Sweep_Angle;
[Wing1_X,Wing1_Y,Wing1_AC,Wing1_S,Wing1_MACloc,Wing1_quarterSweep,Wing1_halfSweep] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Wing1_Y = Wing1_Y - NoseSetbackDist;
Wing1_AC(2) = Wing1_AC(2) - NoseSetbackDist;
Wing1_MACloc = Wing1_MACloc - NoseSetbackDist;

% Tail UAS 1
TailSeperationDist = 0.5199; % [m]
Root_Chord = 0.1448; % [m] 0.160
Tip_Chord = 0.0582; % [m] 0.040
Half_Span = 0.4657/2; % [m] 0.25
tspan1 = 2 * Half_Span;
Sweep_Angle = 16.5827 * (pi/180); Lambdat1 = Sweep_Angle;
[Tail1_X,Tail1_Y,Tail1_AC,Tail1_S,Tail1_MACloc,Tail1_quarterSweep,Tail1_halfSweep] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Tail1_Y = Tail1_Y - TailSeperationDist - NoseSetbackDist;
Tail1_AC(2) = Tail1_AC(2) - TailSeperationDist - NoseSetbackDist;

% % ciurpita.tripod.com/rc/notes/nuetralPt.html
% h0 = 0.25*(Wing1_MACloc(1) - Wing1_MACloc(2)); % Aerodynamic center of wing
% ht = 0.6; % stabilizer efficiency
% Vt = Tail1_S*TailSeperationDist / (Wing1_S*-abs(Wing1_MACloc(1) - Wing1_MACloc(2))); % tail volume
at = 2*pi; % lift curve slope tail
aw = 2*pi; % lift curve slope wing
% NP1 = -h0 + ht*Vt*(at/aw)
% %NP1 = (Tail1_AC(2) - Wing1_AC(2)) / (1 + (Wing1_S / Tail1_S));
% NP1 = NP1 - NoseSetbackDist;

% Spacecraft Design ASU
CMAC = abs(Wing1_MACloc(1) - Wing1_MACloc(2)); CMAC1 = CMAC;
lt = abs(Wing1_AC(2) - Tail1_AC(2)); lt1=lt;
Vt = (lt*Tail1_S) / (Wing1_S*CMAC);
NP1 = Wing1_AC(2) - CMAC*0.6*Vt*(at/aw);


% UAS 2
NoseSetbackDist = 0.100;
width = 0.100; df2 = width;
length = 1.000;
UAS2x = linspace(0,1,100);
t = width/length;
UAS2Fuselage = 5 * t * (0.2969 * sqrt(UAS2x) - 0.1260 * UAS2x - 0.3516 * UAS2x.^2 + 0.2843 * UAS2x.^3 - 0.1015 * UAS2x.^4);
UAS2x = [UAS2x fliplr(UAS2x)] .* length + NoseSetbackDist;
UAS2Fuselage = [UAS2Fuselage -fliplr(UAS2Fuselage)];

% Wing UAS 2
NoseSetbackDist = 0.750; % [m]
Root_Chord = 0.1827; % [m] 0.160
Tip_Chord = 0.1175; % [m] 0.40
Half_Span = 1.3010/2; % [m] 
wspan2 = Half_Span * 2;
Sweep_Angle = -24.6775 * (pi/180); Lambdaw2 = Sweep_Angle;
[Wing2_X,Wing2_Y,Wing2_AC,Wing2_S,Wing2_MACloc,Wing2_quarterSweep,Wing2_halfSweep] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Wing2_Y = Wing2_Y - NoseSetbackDist;
Wing2_AC(2) = Wing2_AC(2) - NoseSetbackDist;
Wing2_MACloc = Wing2_MACloc - NoseSetbackDist;

% Tail UAS 2
TailSeperationDist = -0.300; % [m]
Root_Chord = 0.1054; % [m]
Tip_Chord = 0.0295; % [m]
Half_Span = 0.4847/2; % [m]
tspan2 = Half_Span * 2;
Sweep_Angle = 14.9645 * (pi/180);
[Tail2_X,Tail2_Y,Tail2_AC,Tail2_S,Tail2_MACloc,Tail2_quarterSweep,Tail2_halfSweep] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Tail2_Y = Tail2_Y - TailSeperationDist - NoseSetbackDist;
Tail2_AC(2) = Tail2_AC(2) - TailSeperationDist - NoseSetbackDist;

%NP2 = (Tail2_AC(2) - Wing2_AC(2)) / (1 + (Wing2_S / Tail2_S));
%NP2 = NP2 + Wing2_AC(2);

at = 2*pi; % lift curve slope tail
aw = 2*pi; % lift curve slope wing
% Spacecraft Design ASU
CMAC = abs(Wing2_MACloc(1) - Wing2_MACloc(2)); CMAC2 = CMAC;
lt = abs(Wing2_AC(2) - Tail2_AC(2)); lt2 = lt;
Vt = (lt*Tail2_S) / (Wing2_S*CMAC);
NP2 = Wing2_AC(2) - CMAC*0.6*Vt*(at/aw);

% Combined UAS NP
% Determination of the Nuetral POint for a Box Wing Aircraft [eq 5]

h01 = Wing1_AC(2);
h02 = Wing2_AC(2);
l2 = Wing1_Y(1) - Wing2_Y(1);
c1 = abs(Wing1_MACloc(1) - Wing1_MACloc(2));
c2 = abs(Wing2_MACloc(1) - Wing2_MACloc(2));
CLa1 = 2*pi;
CLa2 = 2*pi;
de_da = 1;
hn = h01 + (l2 + h02*c2 - h01*c1)*(1-de_da)*(CLa2/CLa1) / (1 + (1-de_da)*(CLa2/CLa1));


% UAS 1
figure()
t = tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');
top = tiledlayout(t, 1, 2, 'TileSpacing','compact', 'Padding','compact');

nexttile(top,1)
patch(Wing1_X, Wing1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b'); hold on;
patch(Tail1_X, Tail1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
patch(-Wing1_X, Wing1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
patch(-Tail1_X, Tail1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
scatter([0 0],[Wing1_AC(2) Tail1_AC(2)],20,'bo',"filled");
plot([Wing1_AC(1) Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
plot([-Wing1_AC(1) -Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
plot([Wing1_AC(1) -Wing1_AC(1)], [Wing1_AC(2) Wing1_AC(2)],'k:')
scatter(0,NP1,'kdiamond',"filled");
patch(UAS1Fuselage, -UAS1x, 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b');
grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.1]);
title("UAS 1 : ISR")

% UAS 2
nexttile(top,2)
patch(Wing2_X, Wing2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r'); hold on;
patch(Tail2_X, Tail2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
patch(-Wing2_X, Wing2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
patch(-Tail2_X, Tail2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
scatter([0 0],[Wing2_AC(2) Tail2_AC(2)],20,'ro',"filled");
plot([Wing2_AC(1) Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
plot([-Wing2_AC(1) -Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
plot([Wing2_AC(1) -Wing2_AC(1)], [Wing2_AC(2) Wing2_AC(2)],'k:')
scatter(0,NP2,'kdiamond',"filled");
patch(UAS2Fuselage,-UAS2x, 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'r');
grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.1]);
title("UAS 2 : Payload Drop")

% Combined Vehicle
nexttile(t,2)
patch(Wing1_X, Wing1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b'); hold on;
patch(Tail1_X, Tail1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
patch(Wing2_X, Wing2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
patch(Tail2_X, Tail2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
patch(-Wing1_X, Wing1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
patch(-Tail1_X, Tail1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
patch(-Wing2_X, Wing2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
patch(-Tail2_X, Tail2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
scatter([0 0],[Wing1_AC(2) Tail1_AC(2)],20,'bo',"filled");
scatter([0 0],[Wing2_AC(2) Tail2_AC(2)],20,'ro',"filled");
plot([Wing1_AC(1) Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
plot([Wing2_AC(1) Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
plot([-Wing1_AC(1) -Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
plot([-Wing2_AC(1) -Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
plot([Wing1_AC(1) -Wing1_AC(1)], [Wing1_AC(2) Wing1_AC(2)],'k:')
plot([Wing2_AC(1) -Wing2_AC(1)], [Wing2_AC(2) Wing2_AC(2)],'k:')
%scatter(0,NP1,'kdiamond',"filled");
%scatter(0,NP2,'kdiamond',"filled");
scatter(0,hn,'kdiamond',"filled");

patch(UAS1Fuselage, -UAS1x, 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b');
patch(UAS2Fuselage,-UAS2x, 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'r');
grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.1]);
title("Combined UAS"); hold off;

NP1 = -0.220881;
NP2 = -0.406269;
NPC = -0.336303;


%% longituinal stability %%
% NOTE: xcg/MAC is needed to ensure correct data. h < hn for stability
[figure1,aw,cmt, Cmf, cm0wing, cl_tot1, cd_tot1, alpha1] = UAS1_LongStab(Lambdaw1, Wing1_halfSweep, ...
    wspan1, Wing1_S, df,...
   CMAC1,-Wing1_AC(1,2), -NP1, Tail1_quarterSweep, tspan1, Tail1_S,lt1);

[figure2, aw2, cmc, Cmf2, cm0wing2, cl_tot2, cd_tot2] = UAS2_LongStab(Lambdaw2, Wing2_halfSweep, ...
    wspan2, Wing2_S, df2,...
    CMAC2,-Wing2_AC(1,2), -NP2, Tail2_halfSweep, tspan2, Tail2_S,lt2);

[figure3] = Joined_LongStab(Wing1_S, Wing2_S,wspan1,wspan2, CMAC1, CMAC2,...
    -Wing1_AC(1,2), -Wing2_AC(1,2), aw, aw2, cmt, cmc, cm0wing, cm0wing2, Lambdaw2);

% Lift and drag plots
alpha1 = alpha1.*180/pi;
alpha2 = linspace(-10, 10, 20);

% -----------------------------
% Data (keep your existing vars)
% -----------------------------
cl1_vsp = [-0.4, -0.31, -0.23, -0.15, -0.05, 0.04, 0.13, 0.21, 0.31, 0.4, 0.49, 0.57, 0.65, 0.75, 0.82, 0.9, 1, 1.09, 1.16, 1.25];
cd1_vsp = [0.025, 0.020, 0.016, 0.014, 0.012, 0.0105, 0.012, 0.014, 0.013, 0.016, 0.0185, 0.024, 0.0295, 0.035, 0.042, 0.05, 0.059, 0.069, 0.079, 0.09];  

cl2_vsp = [-0.55, -0.45, -0.35, -0.27, -0.2, -0.1, -0.02, 0.06, 0.15, 0.25, 0.33, 0.40, 0.5, 0.58, 0.65, 0.75, 0.82, 0.9, 1, 1.08];
cd2_vsp = [-0.005, -0.004, -0.003, -0.001, 0, 0.002, 0.004, 0.006, 0.009, 0.012, 0.015, 0.018, 0.02, 0.023, 0.026, 0.03, 0.034, 0.038, 0.042, 0.045];

cl3_vsp = [-0.85, -0.7, -0.55, -0.4, -0.25, -0.1, 0.05, 0.22, 0.38, 0.53, 0.68, 0.81, 0.99, 1.1, 1.25, 1.4, 1.55, 1.7, 1.85, 2];
cd3_vsp = [0.022, 0.02, 0.018, 0.016, 0.015, 0.015, 0.017, 0.019, 0.021, 0.03, 0.04, 0.045, 0.05, 0.06, 0.08, 0.09, 0.105, 0.12, 0.14, .16];

cl_tot3 = cl_tot1 + cl_tot2;
cd_tot3 = cd_tot1 + cd_tot2;


% -----------------------------
% Layout parameters
% -----------------------------
w = 0.35; h = 0.35;
pos1 = [0.08 0.55 w h];
pos2 = [0.57 0.55 w h];
pos3 = [0.32 0.08 w h];

%% =============================
% VEHICLE 1
% =============================
figure('Color','w');

ax1 = axes('Position', pos1);
plot(ax1, alpha2, cl1_vsp);
hold(ax1,'on'); grid(ax1,'on');
plot(ax1, alpha1, cl_tot1);
xlim(ax1, [-10 10]);
xlabel(ax1,'Angle of Attack (deg)');
ylabel(ax1,'C_L');
title(ax1,'Vehicle 1: C_L vs \alpha');
legend(ax1, {'OpenVSP','Calculations'}, 'Location','best');

ax2 = axes('Position', pos2);
hold(ax2,'on'); grid(ax2,'on');
if ~isempty(cd1_vsp)
    plot(ax2, alpha2, cd1_vsp);
end
plot(ax2, alpha1, cd_tot1);
xlabel(ax2,'Angle of Attack (deg)');
ylabel(ax2,'C_D');
title(ax2,'Vehicle 1: C_D vs \alpha');
legend(ax2, {'OpenVSP','Calculations'}, 'Location','best');

ax3 = axes('Position', pos3);
hold(ax3,'on'); grid(ax3,'on');
if ~isempty(cd1_vsp) && ~isempty(cl1_vsp)
    plot(ax3, cd1_vsp, cl1_vsp);
end
plot(ax3, cd_tot1, cl_tot1);
xlabel(ax3,'C_D');
ylabel(ax3,'C_L');
title(ax3,'Vehicle 1: C_L vs C_D');
legend(ax3, {'OpenVSP','Calculations'}, 'Location','best');


%% =============================
% VEHICLE 2
% =============================
figure('Color','w');

ax1 = axes('Position', pos1);
plot(ax1, alpha2, cl2_vsp);
hold(ax1,'on'); grid(ax1,'on');
plot(ax1, alpha1, cl_tot2);
xlim(ax1, [-10 10]);
xlabel(ax1,'Angle of Attack (deg)');
ylabel(ax1,'C_L');
title(ax1,'Vehicle 2: C_L vs \alpha');
legend(ax1, {'OpenVSP','Calculations'}, 'Location','best');

ax2 = axes('Position', pos2);
hold(ax2,'on'); grid(ax2,'on');
if ~isempty(cd2_vsp)
    plot(ax2, alpha2, cd2_vsp);
end
plot(ax2, alpha1, cd_tot2);
xlabel(ax2,'Angle of Attack (deg)');
ylabel(ax2,'C_D');
title(ax2,'Vehicle 2: C_D vs \alpha');
legend(ax2, {'OpenVSP','Calculations'}, 'Location','best');

ax3 = axes('Position', pos3);
hold(ax3,'on'); grid(ax3,'on');
if ~isempty(cd2_vsp) && ~isempty(cl2_vsp)
    plot(ax3, cd2_vsp, cl2_vsp);
end
plot(ax3, cd_tot2, cl_tot2);
xlabel(ax3,'C_D');
ylabel(ax3,'C_L');
title(ax3,'Vehicle 2: C_L vs C_D');
legend(ax3, {'OpenVSP','Calculations'}, 'Location','best');


%% =============================
% VEHICLE 3
% =============================
figure('Color','w');

ax1 = axes('Position', pos1);
plot(ax1, alpha2, cl3_vsp);
hold(ax1,'on'); grid(ax1,'on');
plot(ax1, alpha1, cl_tot3);
xlim(ax1, [-10 10]);
xlabel(ax1,'Angle of Attack (deg)');
ylabel(ax1,'C_L');
title(ax1,'Vehicle 3: C_L vs \alpha');
legend(ax1, {'OpenVSP','Calculations'}, 'Location','best');

ax2 = axes('Position', pos2);
hold(ax2,'on'); grid(ax2,'on');
if ~isempty(cd3_vsp)
    plot(ax2, alpha2, cd3_vsp);
end
plot(ax2, alpha1, cd_tot3);
xlabel(ax2,'Angle of Attack (deg)');
ylabel(ax2,'C_D');
title(ax2,'Vehicle 3: C_D vs \alpha');
legend(ax2, {'OpenVSP','Calculations'}, 'Location','best');

ax3 = axes('Position', pos3);
hold(ax3,'on'); grid(ax3,'on');
if ~isempty(cd3_vsp) && ~isempty(cl3_vsp)
    plot(ax3, cd3_vsp, cl3_vsp);
end
plot(ax3, cd_tot3, cl_tot3);
xlabel(ax3,'C_D');
ylabel(ax3,'C_L');
title(ax3,'Vehicle 3: C_L vs C_D');
legend(ax3, {'OpenVSP','Calculations'}, 'Location','best');