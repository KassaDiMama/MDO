classdef DesignVector
    properties
        % Geometric parameters
        b_half = Const.b_inboard_ref+Const.b_outboard_ref;  
        LE_sweep  = Const.LE_sweep_ref; 
        TR = Const.TR_ref;      
        AR  = Const.AR_ref;      % tip chord [m]

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
            x(1) = obj.b_half;
            x(2) = obj.LE_sweep;
            x(3) = obj.TR;
            x(4) = obj.AR;
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
            obj.b_half = x(1);
            obj.LE_sweep = x(2);
            obj.TR = x(3);
            obj.AR = x(4);
            obj.AU = x(5:10);
            obj.AL = x(11:16);
            obj.Mcr = x(17);
            obj.hcr = x(18);
            
        end
        function str = toString(obj)
            str = sprintf(['DesignVector:\n' ...
                'b_half: %.4f\n' ...
                'LE_sweep: %.4f\n' ...
                'TR: %.4f\n' ...
                'AR: %.4f\n' ...
                'AU: [%s]\n' ...
                'AL: [%s]\n' ...
                'Mcr: %.4f\n' ...
                'hcr: %.2f'], ...
                obj.b_half, obj.LE_sweep, obj.TR, obj.AR, ...
                sprintf('%.4f ', obj.AU), sprintf('%.4f ', obj.AL), ...
                obj.Mcr, obj.hcr);
        end
    end
end