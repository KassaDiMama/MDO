classdef MDA
    properties
        wingDesign WingDesign
        const Const
    end
    
    methods
        
        function obj = MDA(wingDesign, const, W_TO_max, W_ZF)

            obj.wingDesign = wingDesign;
            obj.const = const;
            obj.W_TO_max = W_TO_max;
            obj.W_ZF = W_ZF;
        end
        function [lift_distribution, moment_distribution] = loadsFunc(obj)
                        
            % Wing planform geometry 
            %                x    y     z   chord(m)    twist angle (deg) 
            AC.Wing.Geom = [obj.wingDesign.x_root     obj.wingDesign.y_root     obj.wingDesign.z_root     obj.wingDesign.c_root         obj.wingDesign.twist(1);
                            obj.wingDesign.x_kink     obj.wingDesign.y_kink     obj.wingDesign.z_kink     obj.wingDesign.c_kink         obj.wingDesign.twist(2);
                            obj.wingDesign.x_tip     obj.wingDesign.y_tip     obj.wingDesign.z_tip     obj.wingDesign.c_tip         obj.wingDesign.twist(3)];
            
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
            AC.Aero.V     = obj.wingDesign.V;            % flight speed (m/s)
            AC.Aero.rho   = obj.wingDesign.rho;         % air density  (kg/m3)
            AC.Aero.alt   = obj.wingDesign.hcr;             % flight altitude (m)
            AC.Aero.Re    = obj.wingDesign.Re;        % reynolds number (bqased on mean aerodynamic chord)
            AC.Aero.M     = obj.wingDesign.Mcr;           % flight Mach number 
            AC.Aero.CL    = obj.wingDesign.c_L;          % lift coefficient - comment this line to run the code for given alpha%
            % AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 
            AC.Aero.MaxIterIndex = 150;
            
            tic
            
            % try 
            Res = Q3D_solver(AC);

            y = Res.Wing.Yst;
            cl = Res.Wing.cl;
            cm = Res.Wing.cm_c4;
            chord = Res.Wing.chord;
            L = cl.*chord.*y*0.5*AC.Aero.rho*AC.Aero.V^2;
            M = cm.*chord.*chord.*y*0.5*AC.Aero.rho*AC.Aero.V^2;
            lift_distribution.y = y;
            lift_distribution.L = L;

            moment_distribution.y = y;
            moment_distribution.M = M;

        end
        function [W_TO_max, W_ZF] = structuresFunc(obj, L, M)
            fileName = "optimizing";
            N1 = 0.5;
            N2 = 1;
            %% Create .loading File
            fid = fopen(fileName+'.load', 'w');
        
            if fid == -1
                error('Cannot open loading file for writing.');
            end
            for i = 1:length(y)
                disp(string(y(i)) + " " + string(L(i)) + " " + string(M(i)));
        
                % write the line to file
                fprintf(fid, "%.6f %.6f %.6f\n", y(i), L(i), M(i));
                
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

            
            points_per_side = 46;
            ts = linspace(0,1,points_per_side+1);
            
            
            
            t_upper = flip(ts);     
            t_upper = t_upper(1:end-1);
            y_upper = CST(t_upper, obj.wingDesign.AU);
            
            % Lower surface
            t_lower = ts;
            y_lower = flip(CST(t_lower, obj.wingDesign.AL));
            % Calculate upper surface points using a similar approach
            

            fid = fopen(fileName+'.dat', 'w');

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
                fprintf(fid, string(t)+" "+string(y)+"\n");
            end
            fclose(fid);

            %% Create .init File
            fid = fopen(fileName+'.init', 'w');

            if fid == -1
                error('Cannot open .init file for writing.');
            end
            
            fprintf(fid, string(round(obj.W_TO_max))+" "+string(round(obj.W_ZF))+"\n");
            fprintf(fid, string(obj.const.n_max)+"\n");
            fprintf(fid, "%.2f %.2f %d %d\n", obj.wingDesign.S, obj.wingDesign.b_total, obj.wingDesign.number_of_platforms, obj.wingDesign.number_of_airfoils);
            fprintf(fid, "%.2f %s\n", obj.wingDesign.y_root, fileName);
            fprintf(fid, "%.2f %s\n", obj.wingDesign.y_kink, fileName);
            fprintf(fid, "%.2f %s\n", obj.wingDesign.y_tip, fileName);
            
            fprintf(fid, "%.2f %.2f %.2f %.2f %.2f %.2f\n", ...
            obj.wingDesign.c_root, obj.wingDesign.x_root, obj.wingDesign.y_root, obj.wingDesign.z_root, obj.wingDesign.front_spar_pos,obj.wingDesign.rear_spar_pos);
            fprintf(fid, "%.2f %.2f %.2f %.2f %.2f %.2f\n", ...
            obj.wingDesign.c_kink, obj.wingDesign.x_kink, obj.wingDesign.y_kink, obj.wingDesign.z_kink, obj.wingDesign.front_spar_pos,obj.wingDesign.rear_spar_pos);
            fprintf(fid, "%.2f %.2f %.2f %.2f %.2f %.2f\n", ...
            obj.wingDesign.c_tip, obj.wingDesign.x_tip, obj.wingDesign.y_tip, obj.wingDesign.z_tip, obj.wingDesign.front_spar_pos,obj.wingDesign.rear_spar_pos);
        
            fprintf(fid, "%.2f %.2f\n", obj.wingDesign.start_tank, obj.wingDesign.end_tank);
            fprintf(fid, "%d\n", obj.wingDesign.engine_each_wing);
            fprintf(fid, "%.2f %d\n",obj.wingDesign.engine_location, obj.wingDesign.engine_weight);
        
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            fprintf(fid, "%.5e %.2f %.5e %.5e\n", obj.wingDesign.E_mod, obj.wingDesign.density, obj.wingDesign.tensile_yield, obj.wingDesign.compressive_yield);
            
            fprintf(fid, "%.2f %.2f\n", obj.wingDesign.efficiency_factor, obj.wingDesign.rib_pitch);
            fprintf(fid, "%d\n", 1);
        
            fclose(fid);
            
            
        end
    end
end