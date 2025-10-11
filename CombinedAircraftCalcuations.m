%% Aricraft Nuetral Points
clc; clear; close all;

% UAS 1
width = 0.150;
length = 1.000;
UAS1x = linspace(0,1,100);
t = width/length;
UAS1Fuselage = 5 * t * (0.2969 * sqrt(UAS1x) - 0.1260 * UAS1x - 0.3516 * UAS1x.^2 + 0.2843 * UAS1x.^3 - 0.1015 * UAS1x.^4);
UAS1x = [UAS1x fliplr(UAS1x)] .* length;
UAS1Fuselage = [UAS1Fuselage -fliplr(UAS1Fuselage)];

% Wing UAS 1
NoseSetbackDist = 0.200; % [m]
Root_Chord = 0.240; % [m]
Tip_Chord = 0.060; % [m]
Half_Span = 0.700; % [m]
Sweep_Angle = 30 * (pi/180);
[Wing1_X,Wing1_Y,Wing1_AC,Wing1_S,Wing1_MACloc] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Wing1_Y = Wing1_Y - NoseSetbackDist;
Wing1_AC(2) = Wing1_AC(2) - NoseSetbackDist;
Wing1_MACloc = Wing1_MACloc - NoseSetbackDist;

% Tail UAS 1
TailSeperationDist = 0.550; % [m]
Root_Chord = 0.160; % [m]
Tip_Chord = 0.040; % [m]
Half_Span = 0.250; % [m]
Sweep_Angle = 30 * (pi/180);
[Tail1_X,Tail1_Y,Tail1_AC,Tail1_S,Tail1_MACloc] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Tail1_Y = Tail1_Y - TailSeperationDist - NoseSetbackDist;
Tail1_AC(2) = Tail1_AC(2) - TailSeperationDist - NoseSetbackDist;

% NP1 = -(( Wing1_S * Wing1_AC(2) ) + ( Tail1_S * Tail1_AC(2) )) / ( Wing1_S * Wing1_AC(2) )
%NP1 = (Tail1_AC(2) - Wing1_AC(2)) / (1 + (Wing1_S / Tail1_S));
%NP1 = NP1 + Wing1_AC(2);

% UAS 2
NoseSetbackDist = 0.100;
width = 0.100;
length = 1.000;
UAS2x = linspace(0,1,100);
t = width/length;
UAS2Fuselage = 5 * t * (0.2969 * sqrt(UAS2x) - 0.1260 * UAS2x - 0.3516 * UAS2x.^2 + 0.2843 * UAS2x.^3 - 0.1015 * UAS2x.^4);
UAS2x = [UAS2x fliplr(UAS2x)] .* length + NoseSetbackDist;
UAS2Fuselage = [UAS2Fuselage -fliplr(UAS2Fuselage)];

% Wing UAS 2
NoseSetbackDist = 0.750; % [m]
Root_Chord = 0.160; % [m]
Tip_Chord = 0.040; % [m]
Half_Span = 0.700; % [m]
Sweep_Angle = -10 * (pi/180);
[Wing2_X,Wing2_Y,Wing2_AC,Wing2_S,Wing2_MACloc] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Wing2_Y = Wing2_Y - NoseSetbackDist;
Wing2_AC(2) = Wing2_AC(2) - NoseSetbackDist;
Wing2_MACloc = Wing2_MACloc - NoseSetbackDist;

% Tail UAS 1
TailSeperationDist = -0.500; % [m]
Root_Chord = 0.080; % [m]
Tip_Chord = 0.020; % [m]
Half_Span = 0.200; % [m]
Sweep_Angle = 30 * (pi/180);
[Tail2_X,Tail2_Y,Tail2_AC,Tail2_S,Tail2_MACloc] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle);
Tail2_Y = Tail2_Y - TailSeperationDist - NoseSetbackDist;
Tail2_AC(2) = Tail2_AC(2) - TailSeperationDist - NoseSetbackDist;

%NP2 = (Tail2_AC(2) - Wing2_AC(2)) / (1 + (Wing2_S / Tail2_S));
%NP2 = NP2 + Wing2_AC(2);


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
%scatter(0,NP1,'kdiamond',"filled");
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
%scatter(0,NP2,'kdiamond',"filled");
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

patch(UAS1Fuselage, -UAS1x, 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b');
patch(UAS2Fuselage,-UAS2x, 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'r');
grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.1]);
title("Combined UAS")

