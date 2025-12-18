classdef Optimizer < handle
    properties
        mda MDA
        dvec DesignVector
        wingDesign WingDesign
        x0
        x_normalizer
    end
    
    methods
        function obj = Optimizer(dvec)
            obj.dvec = dvec;
            obj.wingDesign = WingDesign(obj.dvec);
            obj.mda = MDA(obj.wingDesign,Const.W_TO_max_initial,Const.W_ZF_initial);
        end
        function start(obj)
            obj.x0 = obj.dvec.toVector();
            
            objective = @(x) obj.objective_wrapper(x);
            nonlcon = @(x) obj.constraints(x);
            ub = [];
            ub(1) = 52/2 - obj.wingDesign.b_inboard; %b_outboard
            ub(2) = 15;
            ub(3) = 10;
            ub(4) = 3;
            ub(5) = 1;
            ub(6) = 1;
            ub(7) = 1;
            ub(8) = 1;
            ub(9) = 1;
            ub(10) = 1;
            ub(11) = 1;
            ub(12) = 1;
            ub(13) = 1;
            ub(14) = 1;
            ub(15) = 1;
            ub(16) = 1;
            ub(17) = 0.88;
            ub(18) = 13075.92;

            
            lb(1) = 24/2- obj.wingDesign.b_inboard;
            lb(2) = 4;
            lb(3) = 2;
            lb(4) = 1;
            lb(5) = -1;
            lb(6) = -1;
            lb(7) = -1;
            lb(8) = -1;
            lb(9) = -1;
            lb(10) = -1;
            lb(11) = -1;
            lb(12) = -1;
            lb(13) = -1;
            lb(14) = -1;
            lb(15) = -1;
            lb(16) = -1;
            lb(17) = 0.72;
            lb(18) = 10698.48;

            obj.x_normalizer = [];
            obj.x_normalizer(1) = abs(obj.x0(1));
            obj.x_normalizer(2) = abs(obj.x0(2));
            obj.x_normalizer(3) = abs(obj.x0(3));
            obj.x_normalizer(4) = abs(obj.x0(4));
            obj.x_normalizer(5) = 1;
            obj.x_normalizer(6) = 1;
            obj.x_normalizer(7) = 1;
            obj.x_normalizer(8) = 1;
            obj.x_normalizer(9) = 1;
            obj.x_normalizer(10) = 1;
            obj.x_normalizer(11) = 1;
            obj.x_normalizer(12) = 1;
            obj.x_normalizer(13) = 1;
            obj.x_normalizer(14) = 1;
            obj.x_normalizer(15) = 1;
            obj.x_normalizer(16) = 1;
            obj.x_normalizer(17) = abs(obj.x0(17));
            obj.x_normalizer(18) = abs(obj.x0(18));
            x0_normalized = obj.x0./obj.x_normalizer;
            lb_normalized = lb./obj.x_normalizer;
            ub_normalized = ub./obj.x_normalizer;

            A = [];
            b = [];
            Aeq = []; % Equality constraints
            beq = []; % Right-hand side for equality constraints

            % Options (recommended)
            options = optimoptions('fmincon', ...
            'Algorithm','sqp', ...
            'Display','iter-detailed', ...
            'MaxIterations',500, ...
            'MaxFunctionEvaluations',1e6, ...
            'FunctionTolerance',1e-6, ...
            'StepTolerance',1e-6, ...
            'OptimalityTolerance',1e-10, ...
            'DiffMinChange',1e-3, ...
            'DiffMaxChange',1e-2);

            

            [x_opt, fval, exitflag, output] = fmincon(objective, ...
                                          x0_normalized, A, b, Aeq, beq,lb_normalized ,ub_normalized , ...
                                          nonlcon, options);
            disp("Done optimizing!");
        end
        function objective = objective_wrapper(obj,x_normalized)
            x=x_normalized.*obj.x_normalizer;
            obj.dvec = obj.dvec.fromVector(x);
            obj.wingDesign = obj.wingDesign.fromDesignVector(obj.dvec);
            obj.mda.wingDesign = obj.wingDesign;
            % obj.dvec = DesignVector().fromVector(x);
            % obj.wingDesign = WingDesign(obj.dvec);
            % obj.mda.wingDesign = obj.wingDesign;
            objective = -obj.objective_loop();
        end
        function objective = objective_loop(obj)
            lol =obj.mda.MDA_loop(obj.mda.W_TO_max,obj.wingDesign.W_fuel,obj.mda.W_ZF);
            LD_cr = obj.aerodynamicsFunc(obj.mda.W_TO_max,obj.wingDesign.W_fuel);
            eta = obj.performanceFunction();
            objective = obj.objectiveFunc(obj.wingDesign.W_fuel, obj.mda.W_TO_max, LD_cr, eta);
                     

        end

        function [LD_cruise] = aerodynamicsFunc(obj,W_TO_max,W_fuel)
            [CL_wing, CD_wing] = obj.calcCL_CD(W_TO_max,W_fuel);

            q=obj.calculateDesignDynamicPressure();
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
            msg = [
                "Ran objective function that resulted in range = " + string(R)
                "With designVector:"
                obj.dvec.toString()  % column string array
            ];
            
            logMessage(msg, "log.file");
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
                              
            %AC.Wing.eta = [obj.wingDesign.y_root/obj.wingDesign.b_half;obj.wingDesign.y_kink/obj.wingDesign.b_half;obj.wingDesign.y_tip/obj.wingDesign.b_half];  % Spanwise location of the airfoil sections
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
            logMessage([string(datetime('now')) + " | AC details: " + jsonencode(AC)], "log.file");

            Res = Q3D_solver(AC);
            CL_wing = Res.CLwing;
            CD_wing = Res.CDwing;
            if isnan(CD_wing)
                CD_wing = 10;
                disp("CD_wing was NaN")
            end
            logMessage([string(datetime('now')) + " | AC details: " + jsonencode(AC) + " | CL: " + string(CL_wing) + " | CD: " + string(CD_wing)], "log.file");
            % disp("CL_wing = " + string(CL_wing) + ", CD_wing = " + string(CD_wing));
            

        end
        function q = calculateDesignDynamicPressure(obj)
            q=0.5*obj.wingDesign.rho*obj.wingDesign.V^2;
        end
        function [c,ceq] = constraints(obj,x_normalized)
            x=x_normalized.*obj.x_normalizer;
            obj.dvec = obj.dvec.fromVector(x);
            obj.wingDesign = obj.wingDesign.fromDesignVector(obj.dvec);
            obj.mda.wingDesign = obj.wingDesign;

            
            % No inequality constraints
            lol =obj.mda.MDA_loop(obj.mda.W_TO_max,obj.wingDesign.W_fuel,obj.mda.W_ZF);
            c(1) = obj.mda.W_TO_max/obj.wingDesign.S  - Const.W_TO_max_initial/Const.S_ref; % Wing loading constraint  
            ceq = [];
        end
    end
end