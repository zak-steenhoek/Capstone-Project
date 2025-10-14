%% Iterative Code to Solve for Sizing of UAS (MATLAB R2024b)

clc; % Clears Command window
clear; % Clears Workspace
close all; % Closes all open figures

%% Setup Variables for Weight

W_crew = 0; % Weight of crew in lbs
W_payload = 10800; % Weight of payload in lbs

%% Setup Variables for Empty Weight Fraction (Taken from Raymer's Aircraft Design 7th Edition: Table 3.1)

A = 0.93; 
K_vs = 1; % 1.04 if variable sweep, 1.00 for fixed sweep
c = -0.07; % Always negative

%% Weight Fractions At Each Step

% Change amount of steps based on mission profile. 
% For instance, if there is only takeoff, climb, and cruise, 
% you only need to get W3/W0.

% Takeoff and Climb
W_10 = 0.97; % W1/W0
W_21 = 0.985; % W2_W1

% Cruise
R = 9114000; % Range in feet
C_cruise = 0.5/3600; % Specific fuel consumption during cruise from 1/h to 1/s
V = 0.6*(994.8); % Cruise velocity from Mach number to ft/s
L_D_cruise = 16*(0.866); % L/D during cruise

W_32 = exp((-R*C_cruise)/(V*L_D_cruise)); % W3/W2

% Loiter 1
E_1 = 3*3600; % Endurance in loiter 1 from hours to seconds
C_loiter = 0.4/3600; % Specific fuel consumption during loiter 1 from 1/h to 1/s
L_D_loiter = 16; % L/D during loiter 1
W_43 = exp((-E_1*C_loiter)/(L_D_loiter));
W_54 = 0.858;

% Loiter 2
E_2 = (1/3)*3600; % Endurance in loiter 2 from hours to seconds
W_65 = exp((-E_2*C_loiter)/(L_D_loiter));

% Landing
W_76 = 0.995;

% Total weight fraction
W_70 = W_10*W_21*W_32*W_43*W_54*W_65*W_76; % W7/W0
W_f0 = 1.06*(1-W_70); % Wf/W0. Wf/W0 = 0 for battery propulsion systems.

%% Iteration Process for Finding W0

% Setup Variables
W_0_guess = 1000;
error = 0.1;
max_guess = 100;
n = 1;

W_e0 = zeros(max_guess, 1);
W_0 = zeros(max_guess, 1);
W_0(1, 1) = W_0_guess;
guess = zeros(max_guess, 1);

% While loop to calculate W0
while n <= max_guess
    W_e0(n) = A.*K_vs.*((W_0_guess(n)).^c); % We/W0
    W_0(n + 1) = (W_crew + W_payload)./(1 - W_f0 - W_e0(n));
    guess(n, 1) = (W_0(n + 1) - W_0_guess(n))./2;
    if abs(guess(n, 1)) > error
        W_0_guess(n + 1, 1) = W_0(n + 1, 1);
        n = n + 1;
    elseif abs(guess(n, 1)) < error
        fprintf("\nSolution Converged\n\n")
        break
    end
end

% Creates easier to view matrixes in a table
guess_n = (1:1:n)';
guess = nonzeros(guess);
W_e0 = nonzeros(W_e0);
W_e = W_e0.*W_0_guess;
W_0 = nonzeros(W_0);
W_0 = W_0((2:size(W_0, 1)), 1);

Sizing_Table = table(guess_n, W_0_guess, W_e0, W_e, W_0, guess, 'VariableNames', {'Guess', 'W_0, Guess', 'W_e/W_0', 'W_e', 'W_0, Calculated', 'Error'});
% Sizing_Table.("W_0, Guess") = num2str(round(Sizing_Table.("W_0, Guess"), 3));
% Sizing_Table.("W_e/W_0") = num2str(round(Sizing_Table.("W_e/W_0"), 3));
% Sizing_Table.("W_e") = num2str(round(Sizing_Table.("W_e"), 3));
% Sizing_Table.("W_0, Calculated") = num2str(round(Sizing_Table.("W_0, Calculated"), 3));
disp(Sizing_Table)


%% Fsolve for confirmation

W_0_function = @(W_0_f) (W_crew + W_payload)./(1 - W_f0 - (A.*K_vs.*((W_0_f).^c))) - W_0_f;

W_0_f_solve = fsolve(W_0_function, W_0_guess(1, 1)); % Need Optimization Toolbox
fprintf('\nSolution using fsolve: %0.3f lbs \n\n', W_0_f_solve)