classdef WingDesign
    properties
        % Constant Parameters
        b_inboard = 9
        number_of_platforms = 2
        number_of_airfoils = 1
        location_of_airfoil_1 = 0
        
        x_root = 12
        x_kink
        z_kink = 0
        front_spar_pos = 0.1
        rear_spar_pos = 0.9

        start_tank = 0
        end_tank =0.75
        
        engine_each_wing = 1
        engine_location = 0.5
    
        engine_weight = 1900 

        E_mod = 70e3 % in N/mm2 which probably needs to be converted, same for values below
        density = 2800
        tensile_yield = 295
        compressive_yield = 295

        efficiency_factor = 0.96
        rib_pitch =0.5


        % Parameters from DesignVector
        b_outboard = 12.51;     % outboard wing span [m]
        c_root = 8.2;       % root chord [m]
        LE_sweep = 28.55;     % leading-edge sweep [deg]
        c_tip = 1.73;      % tip chord [m]

        AL = ones(1,5) * 0.1;  % lower airfoil shape coefficients
        AU = ones(1,5) * 0.1;  % upper airfoil shape coefficients
        Mcr_ref = 0.8; 
        hcr_ref = 11673.84;    % meters
        Mcr = 0.8;         % initial guess
        hcr = 11673.84;    % initial guess


        % Parameters calculated from above parameters

        S = 10
        b_total
        c_kink = 2
        

        x_tip = 4
        y_tip
        y_kink% todo


        
        

    end

    methods
        %==============================================================
        % Constructor: takes a DesignVector object as input
        %==============================================================
        function obj = WingDesign(dvec)
            arguments
                dvec DesignVector   % input must be a DesignVector object
            end
            % Assign values from DesignVector to WingDesign
            obj.b_outboard = dvec.b_outboard;
            obj.c_root     = dvec.c_root;
            obj.LE_sweep  = dvec.LE_sweep;
            obj.c_tip      = dvec.c_tip;

            obj.AL = dvec.AL;
            obj.AU = dvec.AU;

            obj.Mcr = dvec.Mcr;
            obj.hcr = dvec.hcr;
            obj.x_kink= obj.b_inboard;
            obj.b_total = obj.b_inboard+obj.b_outboard;
            obj.y_tip = obj.b_total;
            obj.y_kink = obj.b_inboard; % todo
        end

        %==============================================================
        function CalculateDesign(obj)
            % TODO: design computations go here
        end

        function calculateSurfaceArea(obj)
            % TODO: geometry calculations go here
        end
    end
end
