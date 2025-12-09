classdef WingDesign < handle
    properties
        % Constant Parameters
        b_inboard = 6.5 % Correct
        number_of_platforms = 2 % Correct
        number_of_airfoils = 3 % Correct
        front_spar_pos = 0.1 % Correct however can be changed later
        rear_spar_pos = 0.9 % Correct however can be changed later
        
        engine_each_wing = 1 % Correct
        engine_location = 0.5 %TBD

        x_root = 0 % Correct
        y_root = 0 % Correct
        z_root = 0 % Correct

        twist = 0; %TBD
        incidence = 0; %TBD
        dihedral = 0; %TBD RADIANS

        E_mod = 70e3 % in N/mm2 which probably needs to be converted, same for values below
        density = 2800 % Correct
        tensile_yield = 295 % Correct but check units
        compressive_yield = 295 % Correct but check units
        
        efficiency_factor = 0.96 %TBD
        rib_pitch =0.5 % Correct


        Mcr_ref = 0.8; %TBD
        hcr_ref = 11673.84; % meters TBD

        engine_weight = 1900 %TBD

        % From Design Vector
        b_outboard     % outboard wing span [m]
        c_root    % root chord [m]
        c_kink
        c_tip     % tip chord [m]

        start_tank
        end_tank
        
        AU % upper airfoil shape coefficients
        AL  % lower airfoil shape coefficients
        Mcr         % initial guess
        hcr    % initial guess
        
        % Calculated
        LE_sweep
        b_total
        S

        x_kink
        x_tip

        y_kink
        y_tip

        z_kink
        z_tip

        rho
        V
        Re

        tank_volume 
        W_fuel
        W_TO_max
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
            obj.c_kink     = dvec.c_kink;
            obj.c_tip      = dvec.c_tip;
            obj.AL         = dvec.AL;
            obj.AU         = dvec.AU;
            obj.Mcr        = dvec.Mcr;
            obj.hcr        = dvec.hcr;


            obj.calculateDesign();
        end

        %==============================================================
        function obj = calculateDesign(obj)
            arguments
                obj
            end

            % Calculate values from current WingDesign
            obj.LE_sweep = obj.calculateLESweep();
            obj.b_total = obj.b_inboard + obj.b_outboard;
            obj.S = obj.calculateSurfaceArea();

            obj.x_kink = obj.c_root - obj.c_kink;
            obj.x_tip = obj.b_total * tan(obj.LE_sweep);

            obj.y_tip = obj.b_total;
            obj.y_kink = obj.b_inboard;

            obj.z_kink = obj.calculateSectionZ(obj.y_kink);
            obj.z_tip = obj.calculateSectionZ(obj.b_total);

            obj.tank_volume = obj.calculateTankVolume();
        end
        
        function LE_sweep = calculateLESweep(obj)
            LE_sweep = atan((obj.c_root-obj.c_kink)/obj.b_inboard);
        end

        function section_z = calculateSectionZ(obj,y)
            section_z = y * tan(obj.dihedral);
        end
        
        function S = calculateSurfaceArea(obj)%, c_r,c_t,LE_sweep, b_outboard)
            % Calculate the surface area of the tapered and kinked wing
            % ckink = c_r - aifoil_X_coord_func(LE_sweep); % chord at kink
            
            S1 = obj.b_inboard/2*(obj.c_root+obj.c_kink);
            S2 = obj.b_outboard/2*(obj.c_kink+obj.c_tip);
        
            S = S1+S2;
        end
        function tank_volume = calculateTankVolume(obj)
            N1 = 0.5;
            N2 = 1;
            function result = CST(t,A)
                cn = t.^N1 .* (1-t).^N2;
                s = 0;
                for i = 0:5
                    s = s + nchoosek(5,i) * t.^i .* (1-t).^(5-i) .* A(i+1);
                end
                result = cn .* s;
            end
            points_per_side = 46;
            ts = linspace(0,1,points_per_side+1);
            
            
            
            t_upper = ts;     
            y_upper = CST(t_upper, obj.AU);
            
            % Lower surface
            t_lower = ts;
            y_lower = -CST(t_lower, obj.AL);

            A_norm = trapz(ts,y_upper) + trapz(ts,y_lower);

            function c = calculateChord(b)
                if b <= obj.b_inboard
                    dc = obj.c_root-obj.c_kink;
                    dc_db = dc/obj.b_inboard;
                    c = obj.c_root - dc_db*b;
                elseif b <= obj.b_total
                    dc = obj.c_kink-obj.c_tip;
                    dc_db = dc/obj.b_outboard;
                    c = obj.c_kink - dc_db*(b-obj.b_inboard);
                else
                    disp("b is higher than b_total")
                end
            end

            bs = linspace(0,obj.b_total,30);
            cs = calculateChord(bs);
            As = A_norm * cs;

            tank_volume=trapz(bs, As);
            disp("Tank Volume"+string(tank_volume));
            % disp(y_lower)

        end
        function dvec = toDesignVector(obj)
            dvec = DesignVector();

            dvec.b_outboard = obj.b_outboard;
            dvec.c_root = obj.c_root;
            dvec.c_kink = obj.c_kink;
            dvec.c_tip = obj.c_tip;
            dvec.AL = obj.AL;
            dvec.AU = obj.AU;
            dvec.Mcr = obj.Mcr;
            dvec.hcr = obj.hcr;

        end
    end
end
