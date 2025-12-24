classdef Initializer < handle
    properties
        % dvec DesignVector
        % wingDesign WingDesign
        % mda MDA
        optimizer Optimizer
        CD_wing_initial
        CL_initial
        q_design_initial
        drag_fus_initial

        W_wing_initial;
        W_ZF_initial
        W_AminusW_initial;

        V_MO_initial 
        V_cr_ref_initial
        wing_tank_volume_initial
        internal_tank_volume
        S_initial
        c_root_initial
        c_kink_initial
        c_tip_initial

        AU
        AL
        AU_lower_bound
        AL_upper_bound
        AU_upper_bound
        AL_lower_bound

        TR_upper_bound

        range_initial
    end
    methods
        function obj = Initializer(dvec)
            arguments
                dvec DesignVector
            end
            % obj.dvec = dvec;
            [obj.AU,obj.AL] = obj.defineCST(Const.airfoil_ref);
            dvec.AU = Const.airfoil_thickness_multiplier*obj.AU;
            dvec.AL = Const.airfoil_thickness_multiplier*obj.AL;
            obj.optimizer = Optimizer(dvec, WingDesign(dvec), obj);
            
            % obj.wingDesign = WingDesign(obj.dvec);
            obj.W_ZF_initial = Const.W_TO_max_initial - Const.W_fuel_cruise_initial;
            % obj.mda = MDA(obj.wingDesign,Const.W_TO_max_initial,obj.W_ZF_initial);
            
            obj.wing_tank_volume_initial = obj.optimizer.wingDesign.calculateWingTankVolume();

            obj.internal_tank_volume = max(0,Const.fuel_weight_max_ref/(Const.rho_fuel)-obj.wing_tank_volume_initial);
            obj.optimizer.wingDesign.internal_tank_volume = obj.internal_tank_volume;
            obj.optimizer.wingDesign.W_fuel = obj.optimizer.wingDesign.calculateFuelWeight();
            [obj.CL_initial, obj.CD_wing_initial] = obj.optimizer.calcCL_CD(Const.W_TO_max_initial,obj.optimizer.wingDesign.W_fuel);
            CDaminusw = obj.CL_initial/Const.LD_initial - obj.CD_wing_initial;
            obj.q_design_initial = obj.optimizer.calculateDesignDynamicPressure();
            obj.drag_fus_initial = CDaminusw * obj.q_design_initial * obj.optimizer.wingDesign.S;
            obj.V_MO_initial = 0.86 * obj.optimizer.wingDesign.a;
            obj.V_cr_ref_initial = obj.optimizer.wingDesign.cruiseSpeed();
            [lift_dist, moment_dist] = obj.optimizer.mda.loadsFunc(Const.W_TO_max_initial,obj.V_MO_initial);
            
            obj.W_wing_initial = obj.optimizer.mda.structuresFunc(lift_dist, moment_dist, Const.W_TO_max_initial, obj.W_ZF_initial);
            
            obj.W_AminusW_initial = obj.W_ZF_initial - obj.W_wing_initial;
            obj.optimizer.mda.W_ZF = obj.W_ZF_initial;
            
            
            
            obj.S_initial = obj.optimizer.wingDesign.S;
            obj.c_root_initial = obj.optimizer.wingDesign.c_root;
            obj.c_kink_initial = obj.optimizer.wingDesign.c_kink;
            obj.c_tip_initial = obj.optimizer.wingDesign.c_tip;
            [obj.AU_lower_bound, obj.AL_upper_bound, obj.AU_upper_bound, obj.AL_lower_bound] = obj.calculateAirfoilBounds();
            
            obj.TR_upper_bound = obj.getTRbounds();

            obj.optimizer.wingDesign.fromDesignVector(obj.optimizer.dvec);
            obj.range_initial = obj.optimizer.objective_loop(); % in meters
            fprintf("Initial range equals: %g km\n",-obj.range_initial/1000);
            
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
                [t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1,N2,AU,AL,"initial_airfoil",flip(ts),CST_order,num_per_side);
                
                % Build one continuous contour (upper forward, lower reversed)
                y = [y_upper; y_lower];
            
                objective = sum((y-airfoil_Y).^2);
            end
            x0 = [AU;AL];
            lb = -50 * ones(1, 2*(CST_order+1));
            ub = -lb;
            
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
        function [AU_lower_bound,AL_upper_bound,AU_upper_bound,AL_lower_bound] =calculateAirfoilBounds(obj)
            function y = CSTcurve(t, A, N1, N2, n)

                % Class function
                C = t.^N1 .* (1 - t).^N2;

                % Shape function 
                S = zeros(size(t));
                for i = 0:n
                    S = S + nchoosek(n, i) .* t.^i .* (1 - t).^(n - i) .* A(i + 1);
                end

                y = C .* S;
            end
            N1 = 0.5;
            N2 = 1;
            CST_order = length(obj.AU) - 1;

            % --- parametric domain ---
            ts = linspace(0, 1, 10000);  % resolution for plotting
            ratios = linspace(0,3,10000);
            % --- compute upper and lower surfaces ---
            for ratio_index = 1:length(ratios)
                ratio = ratios(ratio_index);
                yu = CSTcurve(ts, obj.AU*(1-ratio), N1, N2, CST_order);
                yl = CSTcurve(ts, obj.AL * (1+ratio), N1, N2, CST_order);

                mask_u = yu(2:end-1);
                mask_l = yl(2:end-1);
                res = mask_u > mask_l;
                if sum(res) < size(res,2)
                    % display(ratio)
                    fprintf('Intersection occurred at ratio of %f\n', ratio);
                    AU_lower_bound = (1 - ratios(ratio_index-1));
                    AL_upper_bound = (1 + ratios(ratio_index-1));
                    
                    fprintf('Lower bound: AU = %f, AL = %f\n', AU_lower_bound, AL_upper_bound);
                    
                    fprintf('Lower Ratio(i-1): %f, Upper Ratio (i): %f\n', ratios(ratio_index-1), ratios(ratio_index));
                    fprintf('Difference: %f\n', ratios(ratio_index) - ratios(ratio_index-1));
                    break;
                end
            end
            ratios = linspace(0,20,500);
            for ratio_index = 1:length(ratios)
                ratio = ratios(ratio_index);
                yu = CSTcurve(ts, obj.AU*(1+ratio), N1, N2, CST_order);
                yl = CSTcurve(ts, obj.AL * (1-ratio), N1, N2, CST_order);

                mask_u = yu(2:end-1);
                mask_l = yl(2:end-1);
                res = mask_u > mask_l;
                if max(yu+yl)>0.3
                    % display(ratio)
                    fprintf('Intersection occurred at ratio of %s\n', num2str(ratio));
                    % disp(max(yu+yl));
                    
                    AU_upper_bound = (1 + ratios(ratio_index-1));
                    AL_lower_bound = (1 - ratios(ratio_index-1));
                    break;
                end
            end
            % disp(AU_upper_bound);
            % disp(AL_lower_bound);
            % fprintf("AU: %f",1/max(obj.AU));
            % fprintf("AL: %f",-1/max(obj.AL))
        end
        function [TR_upper_bound] = getTRbounds(obj)
            b_half = Const.b_half_lower_bound;
            b = b_half*2;
            AR = Const.AR_upper_bound;
            LE_sweep = Const.LE_sweep_upper_bound/180*pi;
            b_in = Const.b_inboard_ref;
            TR_upper_bound = (b^2/(AR*b_in*tan(LE_sweep))-b_in)/((b^2/(AR*b_in*tan(LE_sweep)))+b_in+2*(b_half-b_in));
            
        end 
    end
end