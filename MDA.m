classdef MDA
    properties
        wingDesign WingDesign
        const Const
        W_TO_max double
        W_ZF double
    end
    
    methods
        
        function obj = MDA(wingDesign, const, W_TO_max, W_ZF)

            obj.wingDesign = wingDesign;
            obj.const = const;
            obj.W_TO_max = W_TO_max;
            obj.W_ZF = W_ZF;
        end
        function [lift_distribution] = loadsFunc(obj)
                        
            % Wing planform geometry 
            %                x    y     z   chord(m)    twist angle (deg) 
            AC.Wing.Geom = [obj.wingDesign.x_root     obj.wingDesign.y_root     obj.wingDesign.z_root     obj.wingDesign.c_root         obj.wingDesign.twist;
                            obj.wingDesign.x_kink     obj.wingDesign.y_kink     obj.wingDesign.z_kink     obj.wingDesign.c_kink         obj.wingDesign.twist;
                            obj.wingDesign.x_tip     obj.wingDesign.y_tip     obj.wingDesign.z_tip     obj.wingDesign.c_tip         obj.wingDesign.twist];
            
            % Wing incidence angle (degree)
            AC.Wing.inc  = obj.wingDesign.incidence;   
                        
                        
            % Airfoil coefficients input matrix
            %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
            AC.Wing.Airfoils   = [obj.wingDesign.AU; obj.wingDesign.AL;
                                  obj.wingDesign.AU; obj.wingDesign.AL;
                                  obj.wingDesign.AU; obj.wingDesign.AL];
                              
            AC.Wing.eta = [obj.wingDesign.y_root;obj.wingDesign.y_kink;obj.wingDesign.y_tip];  % Spanwise location of the airfoil sections
            
            % Viscous vs inviscid
            AC.Visc  = 1;              % 0 for inviscid and 1 for viscous analysis
            
            % Flight Condition
            AC.Aero.V     = 68;            % flight speed (m/s)
            AC.Aero.rho   = 1.225;         % air density  (kg/m3)
            AC.Aero.alt   = 0;             % flight altitude (m)
            AC.Aero.Re    = 1.14e7;        % reynolds number (bqased on mean aerodynamic chord)
            AC.Aero.M     = 0.2;           % flight Mach number 
            % AC.Aero.CL    = 0.4;          % lift coefficient - comment this line to run the code for given alpha%
            AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 
            AC.Aero.MaxIterIndex = 150;
            
            tic
            
            % try 
            Res = Q3D_solver(AC);
        end
        function [W_TO_max, W_ZF] = structuresFunc(obj)

        end
    end
end