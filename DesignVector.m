classdef DesignVector
    properties
        % Geometric parameters
        b_outboard = 12.51;     % outboard wing span [m]
        c_root     = 8.2;       % root chord [m]
        c_kink     = 6;
        c_tip      = 1.73;      % tip chord [m]

        % Airfoil coefficients Correct whithcomb
        AL = [0.136384890434649	-0.148646719150629	0.056443816239522	-0.365167603195488	-0.120664629058709	-0.613945571520703];  % lower airfoil shape coefficients
        AU = [0.233655457906921	0.079939069830773	0.267462399980753	0.089798155400594	0.277959146200174	0.381599875859377];  % upper airfoil shape coefficients

        % Mach and altitude
        Mcr     = 0.8;         % initial guess
        hcr     = 11673.84;    % initial guess

        % To add
        tank_end

    end
    methods
        function x = toVector(obj)
            x = [];
            x(1) = obj.b_outboard;
            x(2) = obj.c_root;
            x(3) = obj.c_kink;
            x(4) = obj.c_tip;
            x(5) = obj.AU(1);
            x(6) = obj.AU(2);
            x(7) = obj.AU(3);
            x(8) = obj.AU(4);
            x(9) = obj.AU(5);
            x(10) = obj.AU(6);
            x(11) = obj.AL(1);
            x(12) = obj.AL(2);
            x(13) = obj.AL(3);
            x(14) = obj.AL(4);
            x(15) = obj.AL(5);
            x(16) = obj.AL(6);
            x(17) = obj.Mcr;
            x(18) = obj.hcr;
        end
        function obj = fromVector(obj, x)
            obj.b_outboard = x(1);
            obj.c_root = x(2);
            obj.c_kink = x(3);
            obj.c_tip = x(4);
            obj.AU = x(5:9);
            obj.AL = x(10:14);
            obj.Mcr = x(15);
            obj.hcr = x(16);
            
        end
    end
end