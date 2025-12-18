classdef MDA < handle
    properties
        wingDesign WingDesign
        
        W_ZF
        W_TO_max
    end
    
    methods
        
        function obj = MDA(wingDesign,W_TO_max, W_ZF)
            arguments
                wingDesign WingDesign
                W_TO_max
                W_ZF
            end
            obj.wingDesign = wingDesign ;
            obj.W_TO_max = W_TO_max;
            obj.W_ZF = W_ZF;
        end
        function derivative =  W_TO_max_prime(obj)
            delta_W = 1;
            W_TO_plus = obj.W_TO_max+delta_W;
            [lift_distribution_plus, moment_distribution_plus] = obj.loadsFunc(W_TO_plus);
            [W_TO_max_plus, ~] = obj.structuresFunc(lift_distribution_plus, moment_distribution_plus, W_TO_plus, obj.W_ZF);
            derivative = (W_TO_max_plus - obj.W_TO_max) / delta_W;
        end
        
        function new_W = newton(obj)
            [lift_dist, moment_dist] = obj.loadsFunc(obj.W_TO_max);
            [W_TO_calc, W_ZF_calc] = obj.structuresFunc(lift_dist, moment_dist, obj.W_TO_max, obj.W_ZF);
            f_x = obj.W_TO_max - W_TO_calc;
    
           
            d_W_calc_d_W = obj.W_TO_max_prime();
        
            f_prime_x = 1 - d_W_calc_d_W;
        
            if abs(f_prime_x) > 1e-6 % Check for near-zero derivative
                new_W = obj.W_TO_max - (f_x / f_prime_x);
            else

                warning('Derivative is too small. Reverting to simple iteration (W_TO_new = W_TO_calc).');
                new_W = W_TO_calc;
            end
        end

        function obj= MDA_loop(obj,W_TO_max,W_fuel,W_ZF)
            
            W_init = W_TO_max;
            obj.W_TO_max = W_TO_max;
            % obj.W_fuel = W_fuel;
            obj.W_ZF = W_ZF;
            
            [lift_dist, moment_dist] = obj.loadsFunc(W_TO_max);
            obj.W_TO_max, obj.W_ZF = obj.structuresFunc(lift_dist, moment_dist, W_TO_max, W_ZF);
            

            W_TO_old = 0;
            while abs(obj.W_TO_max - W_TO_old) > 4
    
                % Store the current assumed weight before iteration for the convergence check
                W_TO_old = obj.W_TO_max;
                
                % Perform one step of the Newton-Raphson iteration to get a better estimate
                % The 'newton' function already calculates W_TO_calculated for the assumed W_TO_max.
                new_W_TO = obj.newton();
                
                % Update the assumed weight (W_TO_max) with the new, calculated value for the next iteration
                [new_lift_dist, new_moment_dist] = obj.loadsFunc(new_W_TO);
                [new_W_TO, new_ZF] = obj.structuresFunc(new_lift_dist, new_moment_dist, W_TO_max, W_ZF);
                obj.W_TO_max = new_W_TO;
                obj.W_ZF = new_ZF;
                % Optional: display the iteration progress
                % fprintf('Iteration: W_TO_assumed = %.2f, W_TO_new = %.2f, Difference = %.2f\n', W_TO_old, new_W_TO, abs(W_TO_old - new_W_TO));
            end
            [lift_dist, moment_dist] = obj.loadsFunc(obj.W_TO_max);
            [W_TO_final, W_ZF_final] = obj.structuresFunc(lift_dist, moment_dist, obj.W_TO_max, obj.W_ZF);
            obj.W_TO_max = W_TO_final;
            obj.W_ZF = W_ZF_final;
            % fprintf('Finished Iterating: W_TO_initial = %.2f, W_TO_final = %.2f, Difference = %.2f\n', W_init, W_TO_final, abs(W_init - W_TO_final));

        end

        function [lift_distribution, moment_distribution] = loadsFunc(obj, W_TO_max)
                        
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
            AC.Visc  = 0;              % 0 for inviscid and 1 for viscous analysis
            AC.Aero.MaxIterIndex = 150;
            % Flight Condition
            [rho_dont_use,a,T_dont_use] = obj.wingDesign.isa_func();
            Mcritical = Const.V_MO_ref/a;
            Re_corrected = obj.wingDesign.Re/Const.V_MO_ref*obj.wingDesign.V;
            AC.Aero.V     = Const.V_MO_ref;            % flight speed (m/s)
            AC.Aero.rho   = obj.wingDesign.rho;         % air density  (kg/m3)
            AC.Aero.alt   = obj.wingDesign.hcr;             % flight altitude (m)
            AC.Aero.Re    = Re_corrected;        % reynolds number (bqased on mean aerodynamic chord)
            AC.Aero.M     = Mcritical;           % flight Mach number 
            AC.Aero.CL    = obj.wingDesign.calculateCL_critical(W_TO_max);          % lift coefficient - comment this line to run the code for given alpha%
            % AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 

            % tic
            
            % try 
            % disp("Starting Q3D");
            Res = Q3D_solver(AC);
            % disp("Finished Q3D");
            % disp("CL: "+string(AC.Aero.CL))
            % toc
            y     = Res.Wing.Yst(:);          
            cl    = Res.Wing.cl(:);          
            cm    = Res.Wing.cm_c4(:);       
            chord = Res.Wing.chord(:);   
                
            rho = AC.Aero.rho;
            V   = AC.Aero.V;

            lift_distribution.y = y;
            lift_distribution.L = 0.5 * rho * V^2 * chord.*cl;
            moment_distribution.y = y;
            moment_distribution.M = 0.5 * rho * V^2 * chord .* chord .*cm;

            lift_distribution.y = [0;lift_distribution.y; obj.wingDesign.b_total];
            lift_distribution.L = [lift_distribution.L(1);lift_distribution.L;2000];

            moment_distribution.y = [0;moment_distribution.y; obj.wingDesign.b_total];
            moment_distribution.M = [moment_distribution.M(1);moment_distribution.M;-2000];

            L_total = 2 * trapz(lift_distribution.y, lift_distribution.L);
            M_total = 2 * trapz(moment_distribution.y, moment_distribution.M);
                
            % if AC.Visc ==1
            %     disp("Drag Coefficient: "+string(Res.CDwing));
            % end
            % disp("Total lift [N]: "+string(L_total));
            % disp("Cruise lift [N]: "+string(AC.Aero.CL*0.5*rho*V^2*obj.wingDesign.S*2));
            % disp("Cruise lift [kg]: "+string(AC.Aero.CL*0.5*rho*V^2*obj.wingDesign.S*2/9.81));
            % disp("Total moment in Nm"+string(M_total));
        
              
        end
        function [W_TO_max, W_ZF,W_wing] = structuresFunc(obj, lift_distribution, moment_distribution, W_TO_max, W_ZF)
            fileName = "optimizing";
            N1 = 0.5;
            N2 = 1;
            %% Create .loading File
            fid = fopen(fileName+'.load', 'w');
        
            if fid == -1
                error('Cannot open loading file for writing.');
            end
            for i = 1:length(lift_distribution.y)
                % disp(string(lift_distribution.y(i)) + " " + string(lift_distribution.L(i)) + " " + string(moment_distribution.M(i)));
        
                % write the line to file
                if i < length(lift_distribution.y)
                    fprintf(fid, "%.6f %.6f %.6f\n", lift_distribution.y(i)/obj.wingDesign.b_total, lift_distribution.L(i), moment_distribution.M(i));
                else
                    fprintf(fid, "%.6f %.6f %.6f", lift_distribution.y(i)/obj.wingDesign.b_total, lift_distribution.L(i), moment_distribution.M(i));
                end
            end
            fclose(fid);
            %% Create .dat File
            function result = CST(t,A)
                cn = t.^N1 .* (1-t).^N2;
                s = 0;
                for i = 0:5
                    s = s + nchoosek(5,i) * t.^i .* (1-t).^(5-i) .* A(i+1);
                end
                result = cn .* s;
            end

            
            points_per_side = 36;
            ts = linspace(0,1,points_per_side+1);
            
            
            
            t_upper = flip(ts);     
            t_upper = t_upper(1:end-1);
            y_upper = CST(t_upper, obj.wingDesign.AU);
            
            % Lower surface
            t_lower = ts;
            y_lower = CST(t_lower, obj.wingDesign.AL);
            % Calculate upper surface points using a similar approach
            

            fid = fopen('airfoil.dat', 'w');

            if fid == -1
                error('Cannot open .dat file for writing.');
            end
            for index = 1:length(t_upper)
                t = t_upper(index);
                y = y_upper(index);
                fprintf(fid, string(t)+" "+string(y)+"\n");
            end
            for index = 1:length(t_lower)
                t = t_lower(index);
                y = y_lower(index);
                if index < length(t_lower)
                    fprintf(fid, string(t)+" "+string(y)+"\n");
                else
                    fprintf(fid, string(t)+" "+string(y));
                end
            end
            
            fclose(fid);
            
            %% Create .init File
            fid = fopen(fileName+'.init', 'w');

            if fid == -1
                error('Cannot open .init file for writing.');
            end
            
            fprintf(fid, string(round(W_TO_max))+" "+string(round(W_ZF))+"\n");
            fprintf(fid, string(Const.n_max)+"\n");
            fprintf(fid, "%.2f %.2f %d %d\n", obj.wingDesign.S*2, obj.wingDesign.b_total*2, obj.wingDesign.number_of_platforms, obj.wingDesign.number_of_airfoils);
            fprintf(fid, "%.2f %s\n", obj.wingDesign.y_root/obj.wingDesign.b_total, 'airfoil');
            fprintf(fid, "%.2f %s\n", obj.wingDesign.y_tip/obj.wingDesign.b_total, 'airfoil');
            
            fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f\n", ...
            obj.wingDesign.c_root, obj.wingDesign.x_root, obj.wingDesign.y_root, obj.wingDesign.z_root, obj.wingDesign.front_spar_pos,obj.wingDesign.rear_spar_pos);
            fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f\n", ...
            obj.wingDesign.c_kink, obj.wingDesign.x_kink, obj.wingDesign.y_kink, obj.wingDesign.z_kink, obj.wingDesign.front_spar_pos,obj.wingDesign.rear_spar_pos);
            fprintf(fid, "%.4f %.4f %.4f %.4f %.4f %.4f\n", ...
            obj.wingDesign.c_tip, obj.wingDesign.x_tip, obj.wingDesign.y_tip, obj.wingDesign.z_tip, obj.wingDesign.front_spar_pos,obj.wingDesign.rear_spar_pos);
        
            fprintf(fid, "%.4f %.4f\n", obj.wingDesign.start_tank, obj.wingDesign.end_tank);
            fprintf(fid, "%d\n", obj.wingDesign.engine_each_wing);
            fprintf(fid, "%.4f %d\n",obj.wingDesign.engine_location/obj.wingDesign.b_total, obj.wingDesign.engine_weight);
        
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            
            fprintf(fid, "%.2f %.2f\n", obj.wingDesign.efficiency_factor, obj.wingDesign.rib_pitch);
            fprintf(fid, "%d", 0);
        
            fclose(fid);
            
            EMWET optimizing;
            
            fid = fopen(fileName+'.weight','r');  % replace with your file name

            if fid == -1
                error('Cannot open file.');
            end
            
            % Read the first line
            firstLine = fgetl(fid);
            
            % Close the file
            fclose(fid);
            
            % Extract the number from the line using regexp
            massStr = regexp(firstLine, '\d+\.?\d*', 'match');  % match numeric values
            wing_mass = str2double(massStr{1});                 % convert to double
            % disp("wing mass: "+string(wing_mass))
            % disp(obj.wingDesign.LE_sweep*180/pi)
            W_TO_max = Const.W_AminusW+wing_mass+obj.wingDesign.W_fuel;
            W_ZF = Const.W_AminusW+wing_mass;
            W_wing = wing_mass;
        end
    end
end