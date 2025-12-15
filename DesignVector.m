classdef DesignVector
    properties
        % Geometric parameters
        b_outboard = 12.525;   % outboard wing span [m]
        c_root     = 8.2; %11.5 increases weight but also LE sweep       % root chord [m]
        c_kink     = 4.66;      % kink chord [m]
        c_tip      = 1.73;      % tip chord [m]

        % Airfoil coefficients Correct whithcomb
        AL = [-0.2253,-0.1637,-0.0464,-0.4778,0.0741,0.3252];
        AU = [0.2337,0.0800,0.2674,0.0899,0.2779,0.3816];
        % % Airfoil coefficients sc207210
        % AL = [-0.1311,-0.1702,0.0499,-0.4795,0.2237,-0.1444];
        % AU = [0.1603,0.0407,0.2865,-0.1000,0.4355,-0.0540];

        % % Airfoil coefficients 652215
        % AL = [-0.1501,-0.1081,-0.2376,-0.1453,-0.1916,0.0175];
        % AU = [0.1688,0.1787,0.2578,0.2543,0.2285,0.1353];
        
        Mcr     = 0.8;         % initial guess
        hcr     = 11673.84;    % initial guess

        % To add later
        %tank_end

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