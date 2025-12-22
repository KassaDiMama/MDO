classdef Optimizer < handle
    properties
        mda MDA
        dvec DesignVector
        wingDesign WingDesign
        initializer Initializer
        x0
        x_normalizer
    end
    
    methods
        function obj = Optimizer(dvec, wingDesign,initializer)
            obj.dvec =dvec;
            obj.x0 = obj.dvec.toVector();
            obj.wingDesign = wingDesign;
            obj.initializer = initializer;
            obj.mda = MDA(obj.wingDesign,Const.W_TO_max_initial,obj.initializer.W_ZF_initial);
        end
        function start(obj)
            
            
            objective = @(x) obj.objective_wrapper(x);
            nonlcon = @(x) obj.constraints(x);
            ub = [];
            ub(1) = 26/obj.x0(1); %b_outboard
            ub(2) = 40/180*pi/obj.x0(2);
            ub(3) = 0.85/obj.x0(3);
            ub(4) = 12/obj.x0(4);
            ub(5) = obj.initializer.AU_upper_bound;
            ub(6) = obj.initializer.AU_upper_bound;
            ub(7) = obj.initializer.AU_upper_bound;
            ub(8) = obj.initializer.AU_upper_bound;
            ub(9) = obj.initializer.AU_upper_bound;
            ub(10) = obj.initializer.AU_upper_bound;
            ub(11) = obj.initializer.AL_upper_bound;
            ub(12) = obj.initializer.AL_upper_bound;
            ub(13) = obj.initializer.AL_upper_bound;
            ub(14) = obj.initializer.AL_upper_bound;
            ub(15) = obj.initializer.AL_upper_bound;
            ub(16) = obj.initializer.AL_upper_bound;
            ub(17) = 0.88/obj.x0(17);
            ub(18) = 13075.92/obj.x0(18);

            
            lb(1) = 15/obj.x0(1);
            lb(2) = 0/obj.x0(2);
            lb(3) = 0.2/obj.x0(3);
            lb(4) = 6/obj.x0(4);
            lb(5) = obj.initializer.AU_lower_bound;
            lb(6) = obj.initializer.AU_lower_bound;
            lb(7) = obj.initializer.AU_lower_bound;
            lb(8) = obj.initializer.AU_lower_bound;
            lb(9) = obj.initializer.AU_lower_bound;
            lb(10) = obj.initializer.AU_lower_bound;
            lb(11) = obj.initializer.AL_lower_bound;
            lb(12) = obj.initializer.AL_lower_bound;
            lb(13) = obj.initializer.AL_lower_bound;
            lb(14) = obj.initializer.AL_lower_bound;
            lb(15) = obj.initializer.AL_lower_bound;
            lb(16) = obj.initializer.AL_lower_bound;
            lb(17) = 0.72/obj.x0(17);
            lb(18) = 10698.48/obj.x0(18);

   
            x0_normalized = obj.x0./obj.x0;
            

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
            'DiffMinChange',1e-2, ...
            'DiffMaxChange',1e-1);

            disp(string(datetime('now')) + " | Starting optimization...");
            tic;
            
            [x_opt, fval, history, searchdir] = fmincon(objective, ...
                                          x0_normalized, A, b, Aeq, beq,lb ,ub , ...
                                          nonlcon, options);
            toc;
            disp(string(datetime('now')) + " | Finished optimization...");
            disp("Done optimizing!");
            save('fmincon_results.mat', 'x_opt', 'fval', 'history', 'searchdir');
        end
        function objective = objective_wrapper(obj,x_normalized)
            x=x_normalized.*obj.x0;
            obj.dvec = obj.dvec.fromVector(x);
            obj.wingDesign = obj.wingDesign.fromDesignVector(obj.dvec);
            obj.mda.wingDesign = obj.wingDesign;
            % obj.dvec = DesignVector().fromVector(x);
            % obj.wingDesign = WingDesign(obj.dvec);
            % obj.mda.wingDesign = obj.wingDesign;
            range = obj.objective_loop();
            % logmsg("--------------------------------");
            logmsg("At " + string(datetime('now')) + ...
                " calculated range: " + string(range/1000) + " km");
            % logmsg("With wing design:");
            % logmsg(obj.wingDesign.toString());
            % logmsg("--------------------------------");
            objective = -(range/obj.initializer.range_initial);
        end
        function objective = objective_loop(obj)
            lol =obj.mda.MDA_loop(obj.mda.W_TO_max,obj.wingDesign.W_fuel,obj.mda.W_ZF,obj.initializer.W_AminusW_initial,obj.initializer.V_MO_initial);
            LD_cr = obj.aerodynamicsFunc(obj.mda.W_TO_max,obj.wingDesign.W_fuel);
            eta = obj.performanceFunction();
            objective = obj.objectiveFunc(obj.wingDesign.W_fuel, obj.mda.W_TO_max, LD_cr, eta);
                     

        end

        function [LD_cruise] = aerodynamicsFunc(obj,W_TO_max,W_fuel)
            [CL_wing, CD_wing] = obj.calcCL_CD(W_TO_max,W_fuel);

            q=obj.calculateDesignDynamicPressure();
            drag = CD_wing * q *obj.wingDesign.S + obj.initializer.drag_fus_initial/obj.initializer.q_design_initial * q;
            lift = CL_wing * q *obj.wingDesign.S;
            
            LD_cruise = lift/drag;

        end
        function eta = performanceFunction(obj)
            A = (obj.wingDesign.V-obj.initializer.V_cr_ref_initial)^2/(2*70^2);
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
            
            % logMessage(msg, "log.file");
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
            % logMessage([string(datetime('now')) + " | AC details: " + jsonencode(AC)], "log.file");

            Res = Q3D_solver(AC);
            CL_wing = Res.CLwing;
            CD_wing = Res.CDwing;
            if isnan(CD_wing)
                CD_wing = 10;
                disp("CD_wing was NaN")
            end
            % logMessage([string(datetime('now')) + " | AC details: " + jsonencode(AC) + " | CL: " + string(CL_wing) + " | CD: " + string(CD_wing)], "log.file");
            % disp("CL_wing = " + string(CL_wing) + ", CD_wing = " + string(CD_wing));
            

        end
        function q = calculateDesignDynamicPressure(obj)
            q=0.5*obj.wingDesign.rho*obj.wingDesign.V^2;
        end
        function [c,ceq] = constraints(obj,x_normalized)
            x=x_normalized.*obj.x0;
            obj.dvec = obj.dvec.fromVector(x);
            obj.wingDesign = obj.wingDesign.fromDesignVector(obj.dvec);
            obj.mda.wingDesign = obj.wingDesign;

            
            % No inequality constraints
            lol =obj.mda.MDA_loop(obj.mda.W_TO_max,obj.wingDesign.W_fuel,obj.mda.W_ZF,obj.initializer.W_AminusW_initial,obj.initializer.V_MO_initial);
            c(1) = obj.mda.W_TO_max/obj.wingDesign.S  - Const.W_TO_max_initial/obj.initializer.S_initial; % Wing loading constraint  
            ceq = [];
        end
    end
end