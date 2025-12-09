classdef WingDesign
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

        twist = [+4,+2,-5]; %Correct but can be changed later issues
        incidence = 3.2; %degrees Correct
        dihedral = 5; % degrees Correct
        

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
        b_outboard   % outboard wing span [m]
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
        
        rho
        a
        T
        Re   % Reynolds number

        LE_sweep 
        b_total
        S
        MAC

        V
        c_L

        W_fuel
        W_TO_max

        x_kink
        x_tip

        y_kink
        y_tip

        z_kink
        z_tip
        
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
            obj.LE_sweep   = dvec.LE_sweep;
            obj.AL         = dvec.AL;
            obj.AU         = dvec.AU;
            obj.Mcr        = dvec.Mcr;
            obj.hcr        = dvec.hcr;


            obj.calculateDesign(dvec);
        end

        %==============================================================
        function obj = calculateDesign(obj)
            arguments
                obj
            end

            % Calculate values from current WingDesign

            [obj.rho, obj.a, obj.T] = obj.isa_func();

            obj.LE_sweep = obj.calculateLESweep();
            obj.b_total = obj.b_inboard + obj.b_outboard;
            obj.S = obj.calculateSurfaceArea();
            obj.MAC = obj.mac_func();
            obj.V = obj.cruiseSpeed();
            obj.Re = obj.re_func();

            obj.c_L = obj.liftcoef_func();

            obj.x_kink = obj.c_root - obj.c_kink;
            obj.x_tip = obj.b_total * tan(obj.LE_sweep);

            obj.y_tip = obj.b_total;
            obj.y_kink = obj.b_inboard;

            obj.z_kink = obj.calculateSectionZ(obj.y_kink);
            obj.z_tip = obj.calculateSectionZ(obj.b_total);
        end
        
        function LE_sweep = calculateLESweep(obj)
            LE_sweep = arctan((obj.c_root-obj.c_kink)/obj.b_inboard);
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

        function section_z = calculateSectionZ(y)
            section_z = y * tan(deg2rad(obj.dihedral));
        end
        
        function S = calculateSurfaceArea(obj)%, c_r,c_t,LE_sweep, b_outboard)
            % Calculate the surface area of the tapered and kinked wing
            % ckink = c_r - aifoil_X_coord_func(LE_sweep); % chord at kink
            
            S1 = obj.b_inboard/2*(obj.c_root+obj.c_kink);
            S2 = obj.b_outboard/2*(obj.c_kink+obj.c_tip);
        
            S = S1+S2;
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

        function Cl = liftcoef_func(obj)
            %Calculate the required lift coeficient of the aircraft at cruise
            L = sqrt(obj.W_TO_max*(obj.W_TO_max-obj.W_fuel));
            Cl = L/(1/2*obj.rho*obj.V^2*obj.S);
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
