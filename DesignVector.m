classdef DesignVector
    properties
        % Geometric parameters
        b_outboard = 12.51;     % outboard wing span [m]
        c_root     = 8.2;       % root chord [m]
        LE_sweep  = 28.55;     % leading-edge sweep [deg]
        c_tip      = 1.73;      % tip chord [m]

        % Airfoil coefficients (default 5)
        AL = ones(1,5) * 0.1;  % lower airfoil shape coefficients
        AU = ones(1,5) * 0.1;  % upper airfoil shape coefficients

        % Mach and altitude
        Mcr_ref = 0.8; 
        hcr_ref = 11673.84;    % meters
        Mcr     = 0.8;         % initial guess
        hcr     = 11673.84;    % initial guess
    end
    methods
        function CalculateDesign(obj)
            
        end
        function calculateSurfaceArea(obj)

        end
    end
end