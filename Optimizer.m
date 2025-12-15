classdef Optimizer < handle
    properties
        mda MDA
        dvec DesignVector
        wingDesign WingDesign
    end
    
    methods
        function obj = Optimizer(dvec)
            obj.dvec = dvec;
            obj.wingDesign = WingDesign(obj.dvec);
            obj.mda = MDA(obj.wingDesign,Const.W_TO_max_initial,Const.W_ZF_initial);
        end
        function start(obj)
            x0 = obj.dvec.toVector();

            objective = @(x) obj.objective_wrapper(x);


            % Options (recommended)
            options = optimoptions('fmincon',...
                'Algorithm','sqp', ...
                'Display','iter-detailed', ...
                'MaxIterations',500, ...
                'FunctionTolerance',1e-12, ...
                'StepTolerance',1e-12, ...
                'MaxFunctionEvaluations', 1e6, ...
                'OptimalityTolerance',1e-10);
            [x_opt, fval, exitflag, output] = fmincon(objective, ...
                                          x0, A, b, Aeq, beq, lb, ub, ...
                                          nonlcon, options);
        end

        function objective = objective_wrapper(obj,x)
            obj.dvec = DesignVector().fromVector(x);
            obj.wingDesign = WingDesign(obj.dvec);
            obj.mda.wingDesign = obj.wingDesign;
            objective = obj.objective_loop();
        end
        function objective = objective_loop(obj)
            obj.mda.MDA_loop(obj.mda.W_TO_max,obj.wingDesign.W_fuel,obj.mda.W_ZF)
            LD_cr = obj.aerodynamicsFunc(obj.mda.W_TO_max,obj.wingDesign.W_fuel);
            eta = obj.performanceFunction();
            objective = obj.objectiveFunc(obj.wingDesign.W_fuel, obj.mda.W_TO_max, LD_cr, eta);
                     

        end

        function [LD_cruise] = aerodynamicsFunc(obj,W_TO_max,W_fuel)
            [CL_wing, CD_wing] = obj.calcCL_CD(W_TO_max,W_fuel);

            q=0.5*obj.wingDesign.rho*obj.wingDesign.V^2;
            drag = CD_wing * q *obj.wingDesign.S*2 + Const.drag_fus_initial/Const.q_initial * q;
            lift = CL_wing * q *obj.wingDesign.S*2;
            
            LD_cruise = lift/drag;

        end
        function eta = performanceFunction(obj)
            A = (obj.wingDesign.V-Const.Vcr_ref)^2/(2*70^2);
            B = (obj.wingDesign.hcr - Const.hcr_ref)^2/(2*2500^2);
            eta = exp(-A-B);
            
        end
        function objective = objectiveFunc(obj,W_fuel, W_TO_max, LD_cr, eta)
            W_cr_ratio = 1/((W_fuel/W_TO_max-1)/-0.9938);
            CT = Const.CT_bar/eta;
            R = obj.wingDesign.V /CT *LD_cr *log(W_cr_ratio);
            objective = R;
        end
        function [CL_wing,CD_wing] = calcCL_CD(obj,W_TO_max,W_fuel)
            % Wing planform geometry 
            %               x    y     z   chord(m)    twist angle (deg) 
            AC.Wing.Geom = [obj.wingDesign.x_root     obj.wingDesign.y_root     obj.wingDesign.z_root     obj.wingDesign.c_root         obj.wingDesign.twist(1)
                            obj.wingDesign.x_kink     obj.wingDesign.y_kink     obj.wingDesign.z_kink     obj.wingDesign.c_kink         obj.wingDesign.twist(2)
                            obj.wingDesign.x_tip     obj.wingDesign.y_tip     obj.wingDesign.z_tip     obj.wingDesign.c_tip        obj.wingDesign.twist(3)];

            % Wing incidence angle (degree)
            AC.Wing.inc  = obj.wingDesign.incidence;   
                        
                        
            % Airfoil coefficients input matrix
            %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
            AC.Wing.Airfoils   = [obj.wingDesign.AU obj.wingDesign.AL;
                                  obj.wingDesign.AU obj.wingDesign.AL];
                              
            %AC.Wing.eta = [obj.wingDesign.y_root/obj.wingDesign.b_total;obj.wingDesign.y_kink/obj.wingDesign.b_total;obj.wingDesign.y_tip/obj.wingDesign.b_total];  % Spanwise location of the airfoil sections
            AC.Wing.eta = [0;1];
            % Viscous vs inviscid
            AC.Visc  = 1;              % 0 for inviscid and 1 for viscous analysis
            AC.Aero.MaxIterIndex = 150;
            % Flight Condition
            AC.Aero.V     = obj.wingDesign.V;            % flight speed (m/s)
            AC.Aero.rho   = obj.wingDesign.rho;         % air density  (kg/m3)
            AC.Aero.alt   = obj.wingDesign.hcr;             % flight altitude (m)
            AC.Aero.Re    = obj.wingDesign.Re;        % reynolds number (bqased on mean aerodynamic chord)
            AC.Aero.M     = obj.wingDesign.Mcr;           % flight Mach number 
            AC.Aero.CL    = obj.wingDesign.calculateCL_cruise(W_TO_max,W_fuel);          % lift coefficient - comment this line to run the code for given alpha%

            Res = Q3D_solver(AC);

            CL_wing = Res.CLwing;
            CD_wing = Res.CDwing;

        end
    end
end