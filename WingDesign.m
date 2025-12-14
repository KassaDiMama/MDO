classdef WingDesign < handle
    properties
        % Constant Parameters
        b_inboard = 6.5 % Correct
        number_of_platforms = 3 % Correct
        number_of_airfoils = 2 % Correct
        front_spar_pos = 0.2 % Correct however can be changed later
        rear_spar_pos = 0.6 % Correct however can be changed later
        
        engine_each_wing = 1 % Correct
        engine_location

        x_root = 20.5 % Correct CHANGED
        y_root = 0 % Correct
        z_root = 4.1/2 % Correct
        
        start_tank = 0;
        end_tank = 0.85;

        %twist = [+4,+2,-5];%!!!! TURNED -5 at the tip to +2!!! %Correct but can be changed later issues
        twist = [0,0,0]; %CHANGED
        incidence = 3.2; %degrees Correct
        dihedral = 5; % degrees Correct
        

        E_mod = 70e9 % in N/mm2 which probably needs to be converted, same for values below
        density = 2800 % Correct
        tensile_yield = 295e6 % Correct but check units
        compressive_yield = 295e6 % Correct but check units
        
        efficiency_factor = 0.96 % Zed section stiffners, Argued for in report
        rib_pitch = 0.5 % Correct


        Mcr_ref = 0.8; %Correct
        hcr_ref = 11673.84; % meters Corect

        engine_weight = 3705; %kg correct Source: https://en.wikipedia.org/wiki/Rolls-Royce_RB211

        % From Design Vector
        b_outboard   % outboard wing span [m]
        c_root    % root chord [m]
        c_kink
        c_tip     % tip chord [m]

        
        
        AU % upper airfoil shape coefficients
        AL  % lower airfoil shape coefficients
        Mcr         % initial guess
        hcr    % initial guess
        
        % Calculated
        
        rho
        a
        T
        Re   % Reynolds number
        V

        LE_sweep 
        b_total
        S
        MAC



        x_kink
        x_tip

        y_kink
        y_tip

        z_kink
        z_tip

        wing_tank_volume 

        
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

            [obj.rho, obj.a, obj.T] = obj.isa_func();

            
            obj.b_total = obj.b_inboard + obj.b_outboard;
            obj.S = obj.calculateSurfaceArea();
            obj.MAC = obj.mac_func();
            obj.V = obj.cruiseSpeed();
            obj.Re = obj.re_func();


            obj.x_kink = obj.x_root + obj.c_root - obj.c_kink+0.1;
            obj.LE_sweep = obj.calculateLESweep();
            obj.x_tip = obj.x_root + obj.b_total * tan(obj.LE_sweep);

            obj.y_tip = obj.b_total;
            obj.y_kink = obj.b_inboard;

            obj.z_kink = obj.calculateSectionZ(obj.y_kink);
            obj.z_tip = obj.calculateSectionZ(obj.b_total);

            obj.wing_tank_volume = obj.calculateTankVolume();

            obj.engine_location =0.34*obj.b_total; %Correct
        end
        
        function LE_sweep = calculateLESweep(obj)
            LE_sweep = atan((obj.x_kink-obj.x_root)/obj.b_inboard);
        end
        function qc_sweep = calculateQCSweep(obj)
            
            qc_sweep = atan((obj.y_tip-obj.y_root)/(obj.x_tip+0.25*obj.c_tip-obj.x_root-0.25*obj.c_root))*180/pi;
        end
        function MAC = mac_func(obj)
            % Calculate the MAC of the tapered and kinked wing
            tap1 = obj.c_kink/obj.c_root; % taper ratio before kink
            tap2 = obj.c_tip/obj.c_kink; %taper ratio after kink

            mac1 = 2/3 *obj.c_root*(1+tap1+tap1^2)/(1+tap1);
            mac2 = 2/3* obj.c_kink*(1+tap2+tap2^2)/(1+tap2);   
            S1 = 6.5/2*(obj.c_root+obj.c_kink);
            S2 = obj.b_outboard/2*(obj.c_kink+obj.c_tip);

            MAC = (S1*mac1+S2*mac2)/(S1+S2);
        end

        function section_z = calculateSectionZ(obj,y)
            section_z = obj.z_root+ y * tan(deg2rad(obj.dihedral));
        end
        
        function S = calculateSurfaceArea(obj)%, c_r,c_t,LE_sweep, b_outboard)
            % Calculate the surface area of the tapered and kinked wing
            % ckink = c_r - aifoil_X_coord_func(LE_sweep); % chord at kink
            
            S1 = obj.b_inboard/2*(obj.c_root+obj.c_kink);
            S2 = obj.b_outboard/2*(obj.c_kink+obj.c_tip);
        
            S = S1+S2;
        end
        function wing_tank_volume = calculateTankVolume(obj)
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
            ts = linspace(obj.front_spar_pos, obj.rear_spar_pos,points_per_side+1);
            
            
            
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

            wing_tank_volume=trapz(bs, As)*Const.f_tank * 2;
            disp("Tank Volume"+string(wing_tank_volume));

            

            % disp(y_lower)

        end

        function [rho,a, T] = isa_func(obj)
            % ISA_DENSITY_SIMPLE - Calculate ISA air density using piecewise linear approximation
            % Input: h - altitude in meters (0 to 86,000 m)
            % Output: rho - air density in kg/m³
                gamma = 1.4;
                % Convert altitude to kilometers for calculations
                h_km = obj.hcr / 1000;
                
                if h_km < 11
                    % Troposphere (0-11 km)
                    T0 = 288.15;      % Sea level temperature (K)
                    p0 = 101325;      % Sea level pressure (Pa)
                    L = -6.5;         % Temperature lapse rate (K/km)
                    
                    T = T0 + L * h_km;
                    p = p0 * (T / T0).^(-9.80665 / (L / 1000 * 287.05));
                    
                elseif h_km < 20
                    % Lower Stratosphere (11-20 km)
                    T = 216.65;       % Constant temperature (K)
                    p11 = 22632;      % Pressure at 11 km (Pa)
                    
                    p = p11 * exp(-9.80665 * (h_km - 11) * 1000 / (287.05 * T));
                    
                
                    
                else
                    error('Altitude must be between 0 and 32,000 m');
                end
                
                % Calculate density using ideal gas law
                R = 287.05;          % Specific gas constant for air (J/(kg·K))
                rho = p / (R * T); a = sqrt(gamma*R*T);
                
            end
        
        function V = cruiseSpeed(obj)
            V = obj.Mcr * obj.a;
        end

        function Re = re_func(obj)
            % Calculate Reynolds number using the formula Re = (rho * V * h) / mu
                Tref = 273;
                S_ref = 110.4;
                mu0 = 1.716*10^(-5);
                mu = mu0*((obj.T/ Tref)^1.5) * ((Tref+S_ref)/(obj.T+S_ref));
                % Calculate the Reynolds number
                Re = (obj.rho * obj.V * obj.MAC) / mu;
        end

        function Cl = liftcoef_func(obj,W_TO_max, W_fuel)
            %Calculate the required lift coeficient of the aircraft at cruise
            L = sqrt(W_TO_max*(W_TO_max-W_fuel))*9.81;
            % Cl = W_TO_max*9.81/(1/2*obj.rho*obj.V^2*(obj.S*2));
            Cl = L/(1/2*obj.rho*obj.V^2*(obj.S*2));%CHANGED
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
