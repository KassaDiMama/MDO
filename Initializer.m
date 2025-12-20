classdef Initializer < handle
    properties
        dvec DesignVector
        wingDesign WingDesign
        mda MDA
        optimizer Optimizer
        CD_wing_initial
        CL_initial
        q_design_initial
        drag_fus_initial

        W_wing_initial;
        W_ZF_initial
        W_AminusW_initial;

        V_MO_initial 
        
        wing_tank_volume_initial
        internal_tank_volume
        S_initial
        c_root_initial
        c_kink_initial
        c_tip_initial

        AU
        AL

    end
    methods
        function obj = Initializer(dvec)
            arguments
                dvec DesignVector
            end
            obj.dvec = dvec;
            [obj.AU,obj.AL] = obj.defineCST(Const.airfoil_ref);
            obj.dvec.AU = Const.airfoil_thickness_multiplier*obj.AU;
            obj.dvec.AL = Const.airfoil_thickness_multiplier*obj.AL;
            obj.wingDesign = WingDesign(obj.dvec);
            obj.W_ZF_initial = Const.W_TO_max_initial - Const.W_fuel_cruise_initial;
            obj.mda = MDA(obj.wingDesign,Const.W_TO_max_initial,obj.W_ZF_initial);
            obj.optimizer = Optimizer(obj.dvec, obj.wingDesign, obj);
            
            [obj.CL_initial, obj.CD_wing_initial] = obj.optimizer.calcCL_CD(Const.W_TO_max_initial,obj.optimizer.wingDesign.W_fuel);
            CDaminusw = obj.CL_initial/Const.LD_initial - obj.CD_wing_initial;
            obj.q_design_initial = obj.optimizer.calculateDesignDynamicPressure();
            obj.drag_fus_initial = CDaminusw * obj.q_design_initial * obj.wingDesign.S*2;
            obj.V_MO_initial = 0.86 * obj.wingDesign.a;
            
            [lift_dist, moment_dist] = obj.mda.loadsFunc(Const.W_TO_max_initial,obj.V_MO_initial);
            
            obj.W_wing_initial = obj.mda.structuresFunc(lift_dist, moment_dist, Const.W_TO_max_initial, obj.W_ZF_initial);
            
            obj.W_AminusW_initial = obj.W_ZF_initial - obj.W_wing_initial;

            
            obj.wing_tank_volume_initial = obj.wingDesign.calculateWingTankVolume();

            obj.internal_tank_volume = max(0,Const.fuel_weight_max_ref/(Const.rho_fuel)-obj.wing_tank_volume_initial);
            obj.wingDesign.internal_tank_volume = obj.internal_tank_volume;
            obj.S_initial = obj.wingDesign.S;
            obj.c_root_initial = obj.wingDesign.c_root;
            obj.c_kink_initial = obj.wingDesign.c_kink;
            obj.c_tip_initial = obj.wingDesign.c_tip;
            
            
        end
        function [AU, AL] = defineCST(obj,airfoil_name)
            N1 = 0.5;
            N2 = 1;
            CST_order = 5;
            

            % AU = [0.4 0.5 0.5 0.5 0.5 0.1];
            AU = 0.5 * ones(1, CST_order+1);
            % AL = [-0.04 -0.04 -0.04 -0.04 -0.5 -1];
            AL = -0.5 * ones(1, CST_order+1);
            airfoilData = load(airfoil_name);
            airfoil_X = airfoilData(:, 1);
            airfoil_Y = airfoilData(:, 2);

            total_points = length(airfoil_Y);
            num_per_side = (total_points+1)/2;

            function objective = objectiveFunction(x, N1, N2, airfoilData, airfoil_Y,CST_order,num_per_side)
                AU = x(1:CST_order+1);
                AL = x(CST_order+2:2*CST_order+2);
                ts= airfoilData(1:num_per_side, 1);
                [t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1,N2,AU,AL,"test",flip(ts),CST_order,num_per_side);
                
                % Build one continuous contour (upper forward, lower reversed)
                y = [y_upper; y_lower];
            
                objective = sum((y-airfoil_Y).^2);
            end
            x0 = [AU;AL];
            lb = -0.1 * ones(1, 2*(CST_order+1));
            ub = -50*lb;
            
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            nonlcon = [];
            options = optimoptions('fmincon',...
                'Algorithm','sqp', ...
                'Display','iter-detailed', ...
                'MaxIterations',500, ...
                'FunctionTolerance',1e-12, ...
                'StepTolerance',1e-12, ...
                'MaxFunctionEvaluations', 1e6, ...
                'OptimalityTolerance',1e-10);

            objective = @(x) objectiveFunction(x, N1, N2, airfoilData, airfoil_Y,CST_order,num_per_side);
            [x_opt, fval, exitflag, output] = fmincon(objective, ...
                                          x0, A, b, Aeq, beq, lb, ub, ...
                                          nonlcon, options);

            AU = x_opt(1:CST_order+1);
            AL = x_opt(CST_order+2:2*CST_order+2);
        end
    end
end