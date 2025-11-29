%% Tejas Sharma - AEE 478 HW1 Problem 2 Code

clc; % Clears Command window
clear; % Clears Workspace
close all; % Closes all open figures

%%

m = input("\nUAS mass (kg): \nm = ");
z = input("\nAltitude (m): \nz = ");
V_inf = input("\nFlight speed (m/s): \nV_infty = ");
c_root = input("\nChord root length (m): \nc_root = ");
c_tip = input("\nChord tip length (m): \nc_tip = ");
b = input("\nWingspan (m): \nb = ");
Sweep_c4 = input("\nQuarter-chord sweep angle (deg): \nsigma_c/4 = ");

taper_ratio = c_tip./c_root;
c_mgc = (c_tip + c_root)./2;
c_mac = (2/3).*c_root.*((1+taper_ratio+(taper_ratio.^2))./(1+taper_ratio));
Y_bar = (b./6).*((1+2.*taper_ratio)./(1+taper_ratio));
S_guess = (c_root.*(b.*(1+taper_ratio)))./2;

z_star = 8404;
temp_s = 288.15;
density_s = 1.225;
Y = 1.4;
R_air = 287;

Temperature = (1 - ((Y-1)/Y).*(z/z_star)).*temp_s;
mu_air = ((2.791.*(10.^-7)).*(Temperature).^(0.7355));
Density = ((1 - ((Y-1)/Y).*(z/z_star)).^(1/(Y-1))).*density_s;
v_air = mu_air./Density;
Re_c = (V_inf.*c_mac)./(v_air);
a_inf = sqrt(Y.*R_air.*Temperature);
M_inf = V_inf./a_inf;

fprintf("\n\nCurrent Geometric Reference Wing Area # = %.3f m^2", S_guess)
fprintf("\n\nFreestream Mach # = %.3f", M_inf)
fprintf("\n\nReynold's # = %.3f", Re_c)

%%

input("\n\nInput Reynold's # and Mach # into XFOIL for chosen airfoil. \nOnce a text polar accumulation is obtained, hit enter. \n");
PACC_file_name = input("\nInput the absolute path to the PACC text for the airfoil: \n");
Airfoil_Data = readmatrix(string(PACC_file_name));
Airfoil_Data_alpha_zero = find(Airfoil_Data(:, 1) == 0);
Airfoil_Data = cat(1, flip(Airfoil_Data(Airfoil_Data_alpha_zero(4)+1:end, :)), Airfoil_Data(Airfoil_Data_alpha_zero(2):Airfoil_Data_alpha_zero(3)-1, :));
alpha = Airfoil_Data(:, 1);
C_L = Airfoil_Data(:, 2);
C_D = Airfoil_Data(:, 3);
C_Dp = Airfoil_Data(:, 4);
C_M = Airfoil_Data(:, 5);
Top_xtr = Airfoil_Data(:, 6);
Bot_xtr = Airfoil_Data(:, 7);

g = 9.81;
W = m.*g;

%%

Choice = input("\nChoose what to solve for: \n1. Solve for Lift Coefficient and AoA: C_L and AoA \n2. Solve for Reference Wing Area: S \n\n");

if (Choice == 1)
    S = S_guess;
    AR = (b.^2)./S;
    Sweep_LE = rad2deg(atan(tan(deg2rad(Sweep_c4)) + ((1-taper_ratio)/(AR.*(1+taper_ratio)))));
    C_L_AoA = ((2.*W)./(Density.*(V_inf.^2).*S));
    AoA_n = ismembertol(C_L, C_L_AoA, 0.005);
    AoA = alpha(AoA_n);

    fprintf("\nC_L = %.3f at AoA = %.3f degrees \n", C_L_AoA, AoA)

elseif (Choice == 2)
    AoA = input("\nAngle of attack (deg in 0.1 deg intervals): \nalpha = ");

    try
        AoA_n = find(alpha == AoA);
    catch
        error("Input a different angle of attack!!! (0.01 degree intervals)")
    end
    
    C_L_AoA = C_L(AoA_n);
    S = ((2.*W)./(Density.*(V_inf.^2).*C_L_AoA));
    AR = (b.^2)./S;
    Sweep_LE = rad2deg(atan(tan(deg2rad(Sweep_c4)) + ((1-taper_ratio)/(AR.*(1+taper_ratio)))));
    
    fprintf("\n\nS = %.3f m^2 at AoA = %.3f degrees", S, AoA)
    fprintf("\nAR = %.3f \nLeading Edge Sweep = %.3f degrees \n\n", AR, Sweep_LE)

    if (S_guess < S)
        fprintf("Geometric reference wing area is too small. \n - Increase taper ratio (higher c_tip)" + ...
            "\n - Increase wingspan (higher b) \n - Increase flight speed (higher V_inf) \n")
    elseif (S_guess == S)
        fprintf("Geometric reference wing area is good for cruise \n")
    elseif (S_guess > S)
        fprintf("Geometric reference wing area is bigger than needed \n")
    end
end


C_l_max = max(C_L);
C_L_max = 0.9.*C_l_max.*cosd(Sweep_c4);

V_stall = sqrt((2.*W)./(Density.*S.*C_L_max));

C_D_min = min(C_D);
C_l_D_min = C_L(C_D == C_D_min);
if (Sweep_LE < 30)
    e = 1.78.*(1-0.045.*(AR.^(0.68))) - 0.64;
elseif (Sweep_LE > 30 && Sweep_LE < 60)
    e = 4.61.*(1-0.045.*(AR.^(0.68))).*(cosd(Sweep_LE).^(0.15)) - 3.1;
end
K = 1./(pi.*AR.*e); 
C_Di = K.*(C_L_AoA - C_l_D_min).^2;
C_D_AoA = C_D_min + C_Di;

L = (0.5).*Density.*(V_inf.^2).*C_L_AoA.*S;
D = (0.5).*Density.*(V_inf.^2).*C_D_AoA.*S;

fprintf("\nL = %.3f N", L)
fprintf("\nD = %.3f N \n\n", D)


%%

b_x = (linspace(0, b./2, 101))';
c4_line_x = zeros(size(b_x, 1), 1);
c4_line_y = zeros(size(b_x, 1), 1);

c_root_x = (linspace(0, c_root, 101))';
c_root_c4 = c_root_x((c_root_x./c_root) == 0.25);

c_tip_x_start = (b./2).*tand(Sweep_LE);
c_tip_x = (linspace(0, c_tip, 101))';
c_tip_c4 = c_tip_x((c_tip_x./(c_tip)) == 0.25);
c_tip_x = c_tip_x + c_tip_x_start;
c_tip_c4 = c_tip_c4 + c_tip_x_start;

c_mac_x_start = (Y_bar).*tand(Sweep_LE);
c_mac_x = (linspace(0, c_mac, 101))';
c_mac_c4 = c_mac_x((c_mac_x./(c_mac)) == 0.25);
c_mac_x = c_mac_x + c_mac_x_start;
c_mac_c4 = c_mac_c4 + c_mac_x_start;

wing_x_1 = [zeros(size(c_root_x, 1), 1); (-(b./2).*ones(size(c_tip_x, 1), 1)); zeros(size(c_root_x, 1), 1)];
wing_y_1 = [-flip(c_root_x); -c_tip_x; -flip(c_root_x)];

for i = 1:size(b_x, 1)
    c4_line_x(i) = b_x(i);
    c4_line_y(i) = c4_line_x(i).*tand(Sweep_c4) + c_root_c4;
end

%%

figure()
plot(wing_x_1, wing_y_1, "k")
hold on
plot(-wing_x_1, wing_y_1, "k")
hold on
plot(Y_bar.*(ones(size(c_mac_x, 1), 1)), -c_mac_x, "b")
hold on
plot(-Y_bar.*(ones(size(c_mac_x, 1), 1)), -c_mac_x, "b")
hold on
plot(c4_line_x, -c4_line_y, "--r")
hold on
plot(-c4_line_x, -c4_line_y, "--r")
hold on
scatter((b./2), -c_tip_c4, "r", "filled")
hold on
scatter(-(b./2), -c_tip_c4, "r", "filled")
hold on
scatter(Y_bar, -c_mac_c4, "r", "filled")
hold on
scatter(-Y_bar, -c_mac_c4, "r", "filled")
hold on
scatter(0, -c_root_c4, "r", "filled")
xlim([(-(b./2)-0.2) ((b./2)+0.2)])
ylim([(-(b./2)-0.2) ((b./2)+0.2)])

figure()
plot(C_D, C_L, "k")
title("C_{L} vs. C_{D}")
xlabel("Drag coefficicent, C_{D}")
ylabel("Lift coefficicent, C_{L}")

figure()
plot(alpha, C_L, "k")
title("C_{L} vs. {\alpha}")
xlabel("Angle of attack, {\alpha}")
ylabel("Lift coefficicent, C_{L}")

figure()
plot(alpha, C_M, "k")
title("C_{M} vs. {\alpha}")
xlabel("Angle of attack, {\alpha}")
ylabel("Quarter chord moment coefficicent, C_{M}")

figure()
plot(Top_xtr, C_L, "b")
hold on
plot(Bot_xtr, C_L, "r")
title("C_{L} on Top and Bottom Surfaces")
xlabel("X-location along airfoil, xtr")
ylabel("Lift coefficicent, C_{L}")
legend("Top", "Bottom", Location = "best")