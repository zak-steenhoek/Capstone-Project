%% Parameterised UAS Designed
clc; clear; close all;

numberOfItterations = 50000;

Swing1Convergence = 0.1833;
Swing2Convergence = 0.15;

Cr1 = [0.100 0.400];    % [m]
Cr2 = [0.100 0.400];    % [m]
Ct = [0.060 0.200];     % [m]
b = [1.000 1.400];      % [m]
Sweep1 = [20 40];       % [deg]
Sweep2 = [-30 -20];     % [deg] 
x1 = [0.300 0.650];     % [m]
x2 = [-0.300 -0.400];   % [m]

Cr1_t = [0.080 0.250];  % [m]
b1_t = [0.300 0.600];   % [m]
Cr2_t = [0.020 0.200];  % [m]
b2_t = [0.400 0.600];   % [m]
Ct1_t = [0.040 0.060];  % [m]
Sweep1_t = [0 30];      % [deg]
Ct2_t = [0.020 0.040];  % [m]
Sweep2_t = [0 30];      % [deg]

massRatio = [0.830 0.831];

fprintf("Computing %.0f Design Permutations\n\n",numberOfItterations)%length(parameterSpace))
aStart = [mean(Cr1),mean(Cr2),mean(Ct),mean(b),mean(Sweep1),mean(Sweep2),mean(x1),mean(x2),...
    mean(Cr1_t),mean(b1_t),mean(Cr2_t),mean(b2_t),mean(Ct1_t),mean(Sweep1_t),mean(Ct2_t),mean(Sweep2_t),mean(massRatio)];
previousa = aStart;
permutationPercent = 0.05;
previousItteration = 100000;
i = 1;
while i < numberOfItterations
    tic
    a = previousa + ( aStart.*permutationPercent .* ((rand(1,length(aStart)).*2-1)) );
    a = [clip(a(1),Cr1(1),Cr1(2)), clip(a(2),Cr2(1),Cr2(2)), clip(a(3),Ct(1),Ct(2)), clip(a(4),b(1),b(2)), clip(a(5),Sweep1(1),Sweep1(2)), ...
        clip(a(6),Sweep2(1),Sweep2(2)), clip(a(7),x1(1),x1(2)), clip(a(8),x2(1),x2(2)), clip(a(9),Cr1_t(1),Cr1_t(2)), ...
        clip(a(10),b1_t(1),b1_t(2)), clip(a(11),Cr2_t(1),Cr2_t(2)), clip(a(12),b2_t(1),b2_t(2)), clip(a(13),Ct1_t(1),Ct1_t(2)),...
        clip(a(14),Sweep1_t(1),Sweep1_t(2)), clip(a(15),Ct2_t(1),Ct2_t(2)), clip(a(16),Sweep2_t(1),Sweep2_t(2)), clip(a(13),massRatio(1),massRatio(2))];
%     a = [0.248452354142801	0.182663828512092	0.117459967277873	1	26.2489271471051	-24.6774759056559	0.519877481007816	-0.300000000000000	0.144788011777543	0.465664773550571	0.105433108732498	0.484744602166759	0.0581619329085618	16.5827041691186	0.0294809762192670	14.9644789083684	0.830000000000000];

    UAS1_Stability = [0.15 0.35];
    UAS2_Stability = [0.15 0.35];
    Combined_Stability = [0.20 0.30];
    UAS1_mass = 3; % [kg]
    UAS2_mass = UAS1_mass*massRatio; % [kg]
    [NP1, NP2, hn, MAC1, MAC2, Wing1_S, Wing2_S] = CombinedUASInputs(a,UAS1_Stability,UAS2_Stability,Combined_Stability,false);
    
    NP2 = real(NP2); hn = real(hn);

    CG1_for = NP1 + UAS1_Stability(2)*MAC1; CG1_aft = NP1 + UAS1_Stability(1)*MAC1;
    CG2_for = NP2 + UAS2_Stability(2)*MAC2; CG2_aft = NP2 + UAS2_Stability(1)*MAC2;
    CG_Combined_for = hn + Combined_Stability(2)*MAC1; CG_Combined_aft = hn + Combined_Stability(1)*MAC1;
    CG_CombinedCalculated_for = ((CG1_for * UAS1_mass) + (CG2_for * UAS2_mass)) / (UAS1_mass + UAS2_mass);
    CG_CombinedCalculated_aft = ((CG1_aft * UAS1_mass) + (CG2_aft * UAS2_mass)) / (UAS1_mass + UAS2_mass);

    CGDelta(i,1) = min(CG_CombinedCalculated_for, CG_Combined_for);
    CGDelta(i,2) = max(CG_CombinedCalculated_aft, CG_Combined_aft);

    CGParameter = abs(CGDelta(i,1) - CGDelta(i,2));
    S1Parameter = abs(Wing1_S*2 - Swing1Convergence);
    S2Parameter = abs(Wing2_S*2 - Swing2Convergence);
    CostFunction = CGParameter + S1Parameter + S2Parameter;

    if CostFunction < previousItteration
        previousa = a;
        previousItteration = CostFunction;
    end
    ItterationHold(i) = previousItteration;
    
    t(i) = toc;
    fprintf("%.0f/%.0f : %.1f%% Complete : %.2f [ms] : ",i,numberOfItterations,(i/numberOfItterations)*100,t(i)*10^3)
    fprintf("S1: %.2f : S2: %.2f : ",Wing1_S*2,Wing2_S*2)
    if CGDelta(i,1) > CGDelta(i,2)
        fprintf("CG %.1f <-> %.1f [mm]\n",CGDelta(i,1)*1000,CGDelta(i,2)*1000)
    else
        fprintf("Closest Solution : %.1f \n", (abs(CGDelta(i,1) - CGDelta(i,2)))*1000)
    end
    i = i+1;
end
fprintf("\nCalculation Complete : Total Time = %f [s]\n",sum(t))
fprintf("Minimum Cg Difference %f [mm]\n",min(abs(CGDelta(:,1) - CGDelta(:,2)))*1000)

figure()
semilogx(linspace(1,numberOfItterations,numberOfItterations-1),ItterationHold*1000,'k-'); hold on; grid on;
xlabel("Itterations"); ylabel("Cost Function");
title("Parameterization Convergence")

% optimal = find(abs(CGDelta(:,1) - CGDelta(:,2)) == (min(abs(CGDelta(:,1) - CGDelta(:,2)))));
% for i = 1:length(optimal)
%     a = parameterSpace(optimal(i),:);
%
a = previousa;
%%
    % a = [0.248452354142801	0.182663828512092	0.117459967277873	1	26.2489271471051	-24.6774759056559	0.519877481007816	-0.300000000000000	0.144788011777543	0.465664773550571	0.105433108732498	0.484744602166759	0.0581619329085618	16.5827041691186	0.0294809762192670	14.9644789083684	0.830000000000000];
    [NP1, NP2, hn, MAC1, MAC2, Wing1_S, Wing2_S] = CombinedUASInputs(a,UAS1_Stability,UAS2_Stability,Combined_Stability,true);

    NP2 = real(NP2); hn = real(hn);

    CG1_for = NP1 + UAS1_Stability(2)*MAC1; CG1_aft = NP1 + UAS1_Stability(1)*MAC1;
    CG2_for = NP2 + UAS2_Stability(2)*MAC2; CG2_aft = NP2 + UAS2_Stability(1)*MAC2;
    CG_Combined_for = hn + Combined_Stability(2)*MAC1; CG_Combined_aft = hn + Combined_Stability(1)*MAC1;
    CG_CombinedCalculated_for = ((CG1_for * UAS1_mass) + (CG2_for * UAS2_mass)) / (UAS1_mass + UAS2_mass);
    CG_CombinedCalculated_aft = ((CG1_aft * UAS1_mass) + (CG2_aft * UAS2_mass)) / (UAS1_mass + UAS2_mass);

    plot([0 0],[CG1_for CG1_aft],'b-','LineWidth', 4)
    plot([0 0],[CG2_for CG2_aft],'r-','LineWidth', 4)
    plot([0 0],[CG_Combined_for CG_Combined_aft],'k-','LineWidth', 4)
    plot([0 0],[CG_CombinedCalculated_for CG_CombinedCalculated_aft],'m-','LineWidth', 4)
    lgd = legend("","","","","","","","","","","","","","","NP_{UAS Combined}","CG_{UAS 1}","CG_{UAS 2}","CG_{Combined UAS (Combined Masses)}","CG_{Combined UAS (Aerodynamically Determined)}");
    lgd.FontSize = 8;

    set(0, 'Units', 'pixels');
    screenSize = get(0, 'ScreenSize');
    figWidth = 900; figHeight = 700;
    set(gcf, 'Units', 'pixels', 'Position', [(screenSize(3) - figWidth) / 2, (screenSize(4) - figHeight) / 2, figWidth, figHeight]);
% end

fprintf("S1: %.2f : S2: %.2f : ",Wing1_S*2,Wing2_S*2)
    

function [NP1, NP2, hn, CMAC1, CMAC2,Wing1_S,Wing2_S] = CombinedUASInputs(a,UAS1_Stability,UAS2_Stability,Combined_Stability,printresults)
M = 13.90 / 343; % m/s to Mach
Cr1 = a(1);Cr2 = a(2);Ct = a(3);b = a(4);Sweep1 = a(5);Sweep2 = a(6);x1 = a(7);x2 = a(8);Cr1_t = a(9); b1_t = a(10);Cr2_t = a(11);b2_t = a(12);Ct1_t = a(13);Sweep1_t = a(14);Ct2_t = a(15);Sweep2_t = a(16);
% % DEFINE AIRCRAFT % %
% UAS 1
[Wing1_X,Wing1_Y,Wing1_AC,Wing1_S,Wing1_MACloc,Wing1_quarterSweep,Wing1_halfSweep] = MacCode(Cr1, Ct, b/2, Sweep1*(pi/180));
% Tail UAS 1
%Ct1_t = 0.040; % [m]
%Sweep1_t = 30;
[Tail1_X,Tail1_Y,Tail1_AC,Tail1_S,Tail1_MACloc,Tail1_quarterSweep,Tail1_halfSweep] = MacCode(Cr1_t, Ct1_t, b1_t/2, Sweep1_t*(pi/180));
Tail1_Y = Tail1_Y - x1;
Tail1_AC(2) = Tail1_AC(2) - x1;

% UAS 2
[Wing2_X,Wing2_Y,Wing2_AC,Wing2_S,Wing2_MACloc,Wing2_quarterSweep,Wing2_halfSweep] = MacCode(Cr2, Ct, b/2, Sweep2*(pi/180));
NoseSetbackDist = abs(Wing1_Y(2)) + abs(Wing2_Y(2));
Wing2_Y = Wing2_Y - NoseSetbackDist;
Wing2_MACloc = Wing2_MACloc - NoseSetbackDist;
Wing2_AC(2) = Wing2_AC(2) - NoseSetbackDist;
% Tail UAS 2
%Ct2_t = 0.020; % [m]
%Sweep2_t = 30;
[Tail2_X,Tail2_Y,Tail2_AC,Tail2_S,Tail2_MACloc,Tail2_quarterSweep,Tail2_halfSweep] = MacCode(Cr2_t, Ct2_t, b2_t/2, Sweep2_t*(pi/180));
Tail2_Y = Tail2_Y - x2 - NoseSetbackDist;
Tail2_AC(2) = Tail2_AC(2) - x2 - NoseSetbackDist;

Wing1_quarterSweep = Wing1_quarterSweep * pi/180; Tail1_quarterSweep = Tail1_quarterSweep * pi/180;

% % Compute NP % %
% NP1
aw1 = 2*pi; % lift curve slope wing
at1 = 2*pi; % lift curve slope tail
ARwing = b^2 / Wing1_S;
aw1 = (pi*ARwing) * (1 + sqrt(1+ ((1-M^2)*(cos(Wing1_quarterSweep)))*(((pi*ARwing)/(aw1*cos(Wing1_quarterSweep)))^2))).^-1;
ARtail = b1_t^2 / Tail1_S;
at1 = (pi*ARtail) * (1 + sqrt(1+ ((1-M^2)*(cos(Tail1_quarterSweep)))*(((pi*ARtail)/(at1*cos(Tail1_quarterSweep)))^2))).^-1;

CMAC1 = abs(Wing1_MACloc(1) - Wing1_MACloc(2));
lt = abs(Wing1_AC(2) - Tail1_AC(2));
Vt = (lt*Tail1_S) / (Wing1_S*CMAC1);
NP1 = Wing1_AC(2) - CMAC1*0.6*Vt*(at1/aw1);

% NP2
at2 = 2*pi; % lift curve slope tail
aw2 = 2*pi; % lift curve slope wing
ARwing2 = b^2 / Wing2_S;
aw2 = (pi*ARwing2) * (1 + sqrt(1+ ((1-M^2)*(cos(Wing2_quarterSweep)))*(((pi*ARwing2)/(aw2*cos(Wing2_quarterSweep)))^2))).^-1;
ARtail2 = b2_t^2 / Tail2_S;
at2 = (pi*ARtail2) * (1 + sqrt(1+ ((1-M^2)*(cos(Tail2_quarterSweep)))*(((pi*ARtail2)/(at2*cos(Tail2_quarterSweep)))^2))).^-1;

CMAC2 = abs(Wing2_MACloc(1) - Wing2_MACloc(2));
lt = -abs(Wing2_AC(2) - Tail2_AC(2));
Vt = (lt*Tail2_S) / (Wing2_S*CMAC2);
NP2 = Wing2_AC(2) - CMAC2*0.6*Vt*(at2/aw2);

% NP Combined
h01 = Wing1_AC(2); h02 = Wing2_AC(2); l2 = -abs(Wing1_Y(1) - Wing2_Y(1)); de_da = 2*aw1 / (pi*ARwing);
hn = h01 + (l2 + h02*CMAC2 - h01*CMAC1)*(1-de_da)*(aw2/aw1) / (1 + (1-de_da)*(aw2/aw1));


if(printresults)
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
    grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.2]);
    title("UAS 1 : ISR")
    lgd = legend("","","","","1/4 MAC_{W/T}","","","","NP_{UAS 1}");
    lgd.FontSize = 8;

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
    grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.2]);
    title("UAS 2 : Payload Drop")
    lgd = legend("","","","","1/4 MAC_{W/T}","","","","NP_{UAS 2}");
    lgd.FontSize = 8;

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
    plot([Wing1_AC(1) Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
    plot([Wing2_AC(1) Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
    plot([-Wing1_AC(1) -Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
    plot([-Wing2_AC(1) -Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
    plot([Wing1_AC(1) -Wing1_AC(1)], [Wing1_AC(2) Wing1_AC(2)],'k:')
    plot([Wing2_AC(1) -Wing2_AC(1)], [Wing2_AC(2) Wing2_AC(2)],'k:')
    scatter(0,hn,'kdiamond',"filled");
    grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.2]);
    title("Combined UAS"); 
    lgd = legend("","","","","","","","","","","","","","","NP_{UAS Combined}");
    lgd.FontSize = 8;

    set(0, 'Units', 'pixels');
    screenSize = get(0, 'ScreenSize');
    figWidth = 800; figHeight = 700;
    set(gcf, 'Units', 'pixels', 'Position', [(screenSize(3) - figWidth) / 2, (screenSize(4) - figHeight) / 2, figWidth, figHeight]);


    figure()
    patch(Wing1_X, Wing1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b'); hold on;
    patch(Tail1_X, Tail1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
    patch(Wing2_X, Wing2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
    patch(Tail2_X, Tail2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
    patch(-Wing1_X, Wing1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
    patch(-Tail1_X, Tail1_Y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
    patch(-Wing2_X, Wing2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
    patch(-Tail2_X, Tail2_Y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
    plot([Wing1_AC(1) Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
    plot([Wing2_AC(1) Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
    plot([-Wing1_AC(1) -Wing1_AC(1)], [Wing1_MACloc(1) Wing1_MACloc(2)],'k--')
    plot([-Wing2_AC(1) -Wing2_AC(1)], [Wing2_MACloc(1) Wing2_MACloc(2)],'k--')
    plot([Wing1_AC(1) -Wing1_AC(1)], [Wing1_AC(2) Wing1_AC(2)],'k:')
    plot([Wing2_AC(1) -Wing2_AC(1)], [Wing2_AC(2) Wing2_AC(2)],'k:')
    scatter(0,hn,'kdiamond',"filled");
    grid on; axis equal; xlim([-0.8 0.8]); ylim([-1.1 0.2]);
    title("Combined UAS"); 
%     lgd = legend("","","","","","","","","","","","","","","NP_{UAS Combined}");
%     lgd.FontSize = 8;
end
end

function B = clip(A,lower,upper)
B = zeros(1,length(A));
for i = 1:length(A)
    if A(i) < lower
        B(i) = lower;
    elseif A(i) > upper
        B(i) = upper;
    else
        B(i) = A(i);
    end
end

end