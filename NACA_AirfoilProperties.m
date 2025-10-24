function [CD0, CL0, CM0, liftCurveSlope] = NACA_AirfoilProperties(name)
% liftCurveSlope [rad]
% Opperating at Re ~ 200,000
%   MAC = 0.168 [m]
%   V = 20 [m/s]
%   Kinematic Viscosity 1.4207E-5 [m^2/s]

switch name
    case "NACA_0012"
        % http://airfoiltools.com/airfoil/details?airfoil=naca0012h-sa
        CD0 = 0;
        CL0 = 0;
        CM0 = 0;
        liftCurveSlope = 2*pi;
    case "NACA_4412"
        CD0 = 0.01002;
        CL0 = 0.4878;
        CM0 = -0.1078;
        liftCurveSlope = 2*pi;
    otherwise
        CD0 = NaN;
        CL0 = NaN;
        CM0 = NaN;
        liftCurveSlope = NaN;
        disp("Airfoil Does Not Exist")
end

end