
function compareWingDesigns(initial_val, optimizer)
    % Compare Q3D analysis results for initial and optimized wing designs
    % Inputs:
    %   initializer: contains initializer.optimizer.wingDesign
    %   optimizer: contains optimizer.wingDesign
    
    % Extract wing designs
    wingDesign_initial = initial_val.optimizer.wingDesign;
    wingDesign_final = optimizer.wingDesign;
    
    % Calculate results for initial design
    Res_initial = calcQ3D(wingDesign_initial,initial_val.optimizer.mda.W_TO_max,initial_val.optimizer.wingDesign.W_fuel,'Initial design',initial_val);
    
    % Calculate results for final design
    Res_final = calcQ3D(wingDesign_final,optimizer.mda.W_TO_max,optimizer.wingDesign.W_fuel,'final Design',optimizer.initializer);
    
    % Create overlapping plots for drag distribution
    plotOverlappingDrag(Res_initial, Res_final);
    
    % Create overlapping plots for lift distribution
    plotOverlappingLift(Res_initial, Res_final);
end

function Res = calcQ3D(wingDesign,W_TO_max,W_fuel,designName,initializer)
            % Wing planform geometry 
            %               x    y     z   chord(m)    twist angle (deg) 
            AC.Wing.Geom = [wingDesign.x_root     wingDesign.y_root     wingDesign.z_root     wingDesign.c_root         wingDesign.twist(1)
                            wingDesign.x_kink     wingDesign.y_kink     wingDesign.z_kink     wingDesign.c_kink         wingDesign.twist(2)
                            wingDesign.x_tip     wingDesign.y_tip     wingDesign.z_tip     wingDesign.c_tip        wingDesign.twist(3)];

            % Wing incidence angle (degree)
            AC.Wing.inc  = wingDesign.incidence;   
                        
                        
            % Airfoil coefficients input matrix
            %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
            AC.Wing.Airfoils   = [wingDesign.AU wingDesign.AL;
                                  wingDesign.AU wingDesign.AL];
                              
            %AC.Wing.eta = [obj.wingDesign.y_root/obj.wingDesign.b_half;obj.wingDesign.y_kink/obj.wingDesign.b_half;obj.wingDesign.y_tip/obj.wingDesign.b_half];  % Spanwise location of the airfoil sections
            AC.Wing.eta = [0;1];
            % Viscous vs inviscid
            AC.Visc  = 1;              % 0 for inviscid and 1 for viscous analysis
            AC.Aero.MaxIterIndex = 600;
            % Flight Condition
            AC.Aero.V     = wingDesign.V;            % flight speed (m/s)
            AC.Aero.rho   = wingDesign.rho;         % air density  (kg/m3)
            AC.Aero.alt   = wingDesign.hcr;             % flight altitude (m)
            AC.Aero.Re    = wingDesign.Re;        % reynolds number (bqased on mean aerodynamic chord)
            AC.Aero.M     = wingDesign.Mcr;           % flight Mach number 
            AC.Aero.CL    = wingDesign.calculateCL_cruise(W_TO_max,W_fuel);          % lift coefficient - at cruise
            % AC.Aero.CL    = wingDesign.calculateCL_critical(W_TO_max,initializer.V_MO_initial);          % lift coefficient - comment this line to run the code for given alpha%
            % logMessage([string(datetime('now')) + " | AC details: " + jsonencode(AC)], "log.file");
            
            Res = Q3D_solver(AC);
            Res.designName = designName;

            integrand = Res.Wing.cdi .* Res.Wing.chord;
            Res.CDi_total = trapz(Res.Wing.Yst, integrand) / wingDesign.S;
            

            
            
            % logMessage([string(datetime('now')) + " | AC details: " + jsonencode(AC) + " | CL: " + string(CL_wing) + " | CD: " + string(CD_wing)], "log.file");
            % disp("CL_wing = " + string(CL_wing) + ", CD_wing = " + string(CD_wing));
            

end

function plotOverlappingDrag(Res_initial, Res_final)
    % Create overlapped plots of drag distributions for both designs
    % Shows spanwise wing drag coefficients (C∙Cd) with separate curves for:
    % - induced drag contribution
    % - profile+wave drag contribution
    
    % Process initial design
    wingY_initial = Res_initial.Wing.Yst(:);
    wingCdi_initial = Res_initial.Wing.cdi(:);
    secY_initial = Res_initial.Section.Y(:);
    secCd_initial = Res_initial.Section.Cd(:);
    
    % Interpolate section Cd to wing Y locations for initial design
    secCd_on_wing_initial = interp1(secY_initial, secCd_initial, wingY_initial, 'pchip', 'extrap');
    
    % Calculate C*Cd (local chord times local drag coefficient)
    % Get local chord distribution for initial design
    
    chord_initial = Res_initial.Wing.chord(:);
    
    
    % Calculate C*Cd components for initial design
    C_Cdi_initial = chord_initial .* wingCdi_initial;      % Induced drag contribution
    C_Cd_profile_wave_initial = chord_initial .* secCd_on_wing_initial;  % Profile+wave drag contribution
    C_Cd_total_initial = C_Cdi_initial + C_Cd_profile_wave_initial;      % Total
    
    % Process final design
    wingY_final = Res_final.Wing.Yst(:);
    wingCdi_final = Res_final.Wing.cdi(:);
    secY_final = Res_final.Section.Y(:);
    secCd_final = Res_final.Section.Cd(:);
    
    % Interpolate section Cd to wing Y locations for final design
    secCd_on_wing_final = interp1(secY_final, secCd_final, wingY_final, 'pchip', 'extrap');
    
    % Get local chord distribution for final design
    
    chord_final = Res_final.Wing.chord(:);
    
    
    % Calculate C*Cd components for final design
    C_Cdi_final = chord_final .* wingCdi_final;          % Induced drag contribution
    C_Cd_profile_wave_final = chord_final .* secCd_on_wing_final;  % Profile+wave drag contribution
    C_Cd_total_final = C_Cdi_final + C_Cd_profile_wave_final;      % Total
    
    % Create figure for drag coefficients
    figure('Name', 'Spanwise Wing Drag Coefficients (C∙Cd) @ Design Point', 'Position', [100, 100, 1400, 600]);
    
    % ===== SUBPLOT 2: Overlay of Induced Drag Contribution =====
    subplot(1, 2, 1);  % Changed from 2,3,2 to 1,2,1
    hold on;
    
    % Initial design - Induced drag
    area(wingY_initial, C_Cdi_initial, 'FaceColor', 'b', 'FaceAlpha', 0.3, ...
        'EdgeColor', 'b', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Induced - %s', Res_initial.designName));
    
    % Final design - Induced drag
    area(wingY_final, C_Cdi_final, 'FaceColor', 'r', 'FaceAlpha', 0.3, ...
        'EdgeColor', 'r', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Induced - %s', Res_final.designName));
    
    % Add line plots on top of areas
    plot(wingY_initial, C_Cdi_initial, 'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(wingY_final, C_Cdi_final, 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    xlabel('Spanwise Location (m)');
    ylabel('C∙Cd_i (m)');
    title('Induced Drag Contribution (C∙Cd_i)');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    box on;
    hold off;
    
    % ===== SUBPLOT 3: Overlay of Profile+Wave Drag Contribution =====
    subplot(1, 2, 2);  % Changed from 2,3,3 to 1,2,2
    hold on;
    
    % Initial design - Profile+Wave drag
    area(wingY_initial, C_Cd_profile_wave_initial, 'FaceColor', 'b', 'FaceAlpha', 0.3, ...
        'EdgeColor', 'b', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Profile+Wave - %s', Res_initial.designName));
    
    % Final design - Profile+Wave drag
    area(wingY_final, C_Cd_profile_wave_final, 'FaceColor', 'r', 'FaceAlpha', 0.3, ...
        'EdgeColor', 'r', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Profile+Wave - %s', Res_final.designName));
    
    % Add line plots on top of areas
    plot(wingY_initial, C_Cd_profile_wave_initial, 'b--', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(wingY_final, C_Cd_profile_wave_final, 'r--', 'LineWidth', 2, 'HandleVisibility', 'off');
    
    xlabel('Spanwise Location (m)');
    ylabel('C∙Cd_{p+w} (m)');
    title('Profile+Wave Drag Contribution (C∙Cd_{p+w})');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    box on;
    hold off;
    
    % ===== SUBPLOT 6: Drag Ratio Improvement (COMMENTED OUT) =====
    %{
    subplot(2, 3, 6);
    
    % Calculate percentage improvements
    improvement_C_Cdi = 100 * (integrated_C_Cdi_initial - integrated_C_Cdi_final) / integrated_C_Cdi_initial;
    improvement_C_Cd_pw = 100 * (integrated_C_Cd_pw_initial - integrated_C_Cd_pw_final) / integrated_C_Cd_pw_initial;
    improvement_total = 100 * (integrated_total_initial - integrated_total_final) / integrated_total_initial;
    
    improvements = [improvement_C_Cdi, improvement_C_Cd_pw, improvement_total];
    categories_improvement = {'Induced Drag', 'Profile+Wave', 'Total'};
    
    bar(improvements, 'FaceColor', [0.2, 0.6, 0.2], 'EdgeColor', 'k');
    
    % Add value labels
    for i = 1:length(improvements)
        if improvements(i) > 0
            text_color = 'g';
            prefix = '▼';
        else
            text_color = 'r';
            prefix = '▲';
        end
        text(i, improvements(i), sprintf('%s%.1f%%', prefix, abs(improvements(i))), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 9, 'FontWeight', 'bold', 'Color', text_color);
    end
    
    set(gca, 'XTick', 1:length(categories_improvement), 'XTickLabel', categories_improvement);
    ylabel('Improvement (%)');
    title('Drag Reduction Percentage');
    grid on;
    box on;
    yline(0, 'k--', 'LineWidth', 1);  % Zero line reference
    %}
    
    % Add overall title
    sgtitle('Spanwise Wing Drag Coefficients (C∙Cd) @ Design Point: Initial vs Optimized Design', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Display summary in command window
    fprintf('\n=== Drag Analysis Summary ===\n');
    
    % Calculate integrated drag values for summary
    integrated_C_Cdi_initial = trapz(wingY_initial, C_Cdi_initial);
    integrated_C_Cd_pw_initial = trapz(wingY_initial, C_Cd_profile_wave_initial);
    integrated_total_initial = integrated_C_Cdi_initial + integrated_C_Cd_pw_initial;
    
    integrated_C_Cdi_final = trapz(wingY_final, C_Cdi_final);
    integrated_C_Cd_pw_final = trapz(wingY_final, C_Cd_profile_wave_final);
    integrated_total_final = integrated_C_Cdi_final + integrated_C_Cd_pw_final;
    
    % Calculate improvements
    improvement_C_Cdi = 100 * (integrated_C_Cdi_initial - integrated_C_Cdi_final) / integrated_C_Cdi_initial;
    improvement_C_Cd_pw = 100 * (integrated_C_Cd_pw_initial - integrated_C_Cd_pw_final) / integrated_C_Cd_pw_initial;
    improvement_total = 100 * (integrated_total_initial - integrated_total_final) / integrated_total_initial;
    
    fprintf('Integrated C∙Cd values:\n');
    fprintf('  %s: Induced = %.4f m², Profile+Wave = %.4f m², Total = %.4f m²\n', ...
        Res_initial.designName, integrated_C_Cdi_initial, integrated_C_Cd_pw_initial, integrated_total_initial);
    fprintf('  %s: Induced = %.4f m², Profile+Wave = %.4f m², Total = %.4f m²\n', ...
        Res_final.designName, integrated_C_Cdi_final, integrated_C_Cd_pw_final, integrated_total_final);
    fprintf('\nImprovements:\n');
    fprintf('  Induced Drag: %.1f%% reduction\n', improvement_C_Cdi);
    fprintf('  Profile+Wave Drag: %.1f%% reduction\n', improvement_C_Cd_pw);
    fprintf('  Total Drag: %.1f%% reduction\n', improvement_total);
    fprintf('=============================\n\n');
end
function plotOverlappingLift(Res_initial, Res_final)
    % Create overlapped plots of lift distributions for both designs
    
    figure('Name', 'Lift Distribution Comparison at design condition', 'Position', [100, 100, 1200, 600]);
    
    % Subplot 1: Lift distribution comparison
   
    hold on;
    
    % Plot initial design
    
    plot(Res_initial.Wing.Yst, Res_initial.Wing.ccl, 'b-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Lift Coef - %s', Res_initial.designName));
    
    
    % Plot final design
    
    plot(Res_final.Wing.Yst, Res_final.Wing.ccl, 'r-', 'LineWidth', 2, ...
        'DisplayName', sprintf('Lift Coef - %s', Res_final.designName));
    
    
    
    xlabel('Spanwise Location');
    ylabel('Lift Coefficient');
    title('Lift Distribution Comparison');
    legend show;
    grid on;
    axis equal;
    
    % Get current axis limits to see what we have
    x_limits = xlim;
    y_limits = ylim;
    
    % Determine which axis needs more padding
    x_range = x_limits(2) - x_limits(1);
    y_range = y_limits(2) - y_limits(1);
    
    % Add 10% padding to both axes
    padding_factor = 0.1;
    x_padding = x_range * padding_factor;
    y_padding = y_range * padding_factor;
    
    % Set new limits with padding
    xlim([x_limits(1), x_limits(2)]);
    ylim([y_limits(1) - y_padding, y_limits(2) + y_padding]);
    
    hold off;
    
    
    
end
% clear all
% close all
% clc
% dvec1 = DesignVector();
% initial_values = Initializer(dvec1);

% clear all
%load files
initializer = load("fmincon_2026-01-07_11-55-38\initializer2026-01-07_11-51-48.mat").initializer;
final_x_normalized = load("fmincon_2026-01-07_11-55-38\final.mat","x").x;

% initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;
% final_x_normalized = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\final.mat","x").x;
%

optimizer= initializer.optimizer;

final_x = final_x_normalized.*optimizer.x0;
dvec = DesignVector().fromVector(final_x);
wingDesign = WingDesign(dvec);
optimizer.dvec = dvec;
optimizer.wingDesign = wingDesign;
optimizer.mda.wingDesign = wingDesign;

optimizer.mda.MDA_loop(Const.W_TO_max_initial,Const.W_fuel_cruise_initial,initializer.W_ZF_initial,initializer.W_AminusW_initial,initializer.V_MO_initial);
W_fuel = optimizer.mda.W_TO_max-optimizer.mda.W_ZF;
[CL_wing, CD_wing]=optimizer.calcCL_CD(optimizer.mda.W_TO_max,W_fuel); 
display(CL_wing)
display(CD_wing)
LD = optimizer.aerodynamicsFunc(optimizer.mda.W_TO_max,W_fuel);
display(LD)
            % Wing planform geometry 
