
function compareWingDesigns(initial_val,initializer, optimizer)
    % Compare Q3D analysis results for initial and optimized wing designs
    % Inputs:
    %   initializer: contains initializer.optimizer.wingDesign
    %   optimizer: contains optimizer.wingDesign
    
    % Extract wing designs
    wingDesign_initial = initial_val.optimizer.wingDesign;
    wingDesign_final = optimizer.wingDesign;
    
    % Calculate results for initial design
    Res_initial = runQ3DAnalysis(wingDesign_initial, initializer, 'Initial Design',initial_val.optimizer.mda.W_TO_max)
    
    % Calculate results for final design
    Res_final = runQ3DAnalysis(wingDesign_final, initializer, 'Final Design',optimizer.mda.W_TO_max)
    
    % Create overlapping plots for drag distribution
    % plotOverlappingDrag(Res_initial, Res_final);
    
    % Create overlapping plots for lift distribution
    plotOverlappingLift(Res_initial, Res_final);
end

function Res = runQ3DAnalysis(wingDesign, initializer, designName, W_TO_max)
    % Run Q3D analysis for a given wing design
    fprintf('Running Q3D analysis for %s...\n', designName);

    
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
                      
    AC.Wing.eta = [0;1];
    
    % Viscous vs inviscid
    AC.Visc  = 1;              % 0 for inviscid and 1 for viscous analysis
    AC.Aero.MaxIterIndex = 150;
    
    % Flight Condition
    [rho_dont_use, a, T_dont_use] = wingDesign.isa_func();
    Mcritical = initializer.V_MO_initial / a;
    Re_corrected = wingDesign.Re / initializer.V_MO_initial * wingDesign.V;
    
    AC.Aero.V     = initializer.V_MO_initial;            % flight speed (m/s)
    AC.Aero.rho   = wingDesign.rho;         % air density  (kg/m3)
    AC.Aero.alt   = wingDesign.hcr;             % flight altitude (m)
    AC.Aero.Re    = Re_corrected;        % reynolds number (based on mean aerodynamic chord)
    AC.Aero.M     = Mcritical;           % flight Mach number 
    
    % Note: You need to define 'const.W_TO_max_initial' or pass it as parameter
    % For now, I'll assume it's available in the workspace or we'll calculate CL differently
    % You may need to adjust this line based on your actual Const structure:
    AC.Aero.CL = wingDesign.calculateCL_cruise(W_TO_max, initializer.V_MO_initial);
    
    % Run Q3D solver
    Res = Q3D_solver(AC);
    
    % Store design name for plotting
    Res.designName = designName;
end

function plotOverlappingDrag(Res_initial, Res_final)
    % Create overlapped plots of drag distributions for both designs
    
    % Process initial design
    wingY_initial = Res_initial.Wing.Yst(:);
    wingCdi_initial = Res_initial.Wing.cdi(:);
    secY_initial = Res_initial.Section.Y(:);
    secCd_initial = Res_initial.Section.Cd(:);
    
    % Interpolate section Cd to wing Y locations for initial design
    secCd_on_wing_initial = interp1(secY_initial, secCd_initial, wingY_initial, 'pchip', 'extrap');
    CD_total_initial = wingCdi_initial + secCd_on_wing_initial;
    
    % Process final design
    wingY_final = Res_final.Wing.Yst(:);
    wingCdi_final = Res_final.Wing.cdi(:);
    secY_final = Res_final.Section.Y(:);
    secCd_final = Res_final.Section.Cd(:);
    
    % Interpolate section Cd to wing Y locations for final design
    secCd_on_wing_final = interp1(secY_final, secCd_final, wingY_final, 'pchip', 'extrap');
    CD_total_final = wingCdi_final + secCd_on_wing_final;
    
    % Create figure for drag coefficients
    figure('Name', 'Drag Coefficients Comparison', 'Position', [100, 100, 1200, 800]);
    
    % Subplot 1: All components for initial design
    subplot(2, 2, 1);
    hold on;
    plot(Res_initial.Wing.Yst, Res_initial.Wing.cdi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Wing CDi');
    plot(Res_initial.Section.Y, Res_initial.Section.Cd, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Section CD');
    plot(wingY_initial, CD_total_initial, 'g-', 'LineWidth', 2, 'DisplayName', 'Total CD');
    xlabel('Spanwise Location');
    ylabel('Drag Coefficient');
    title(sprintf('Drag Components - %s', Res_initial.designName));
    legend show;
    grid on;
    hold off;
    
    % Subplot 2: All components for final design
    subplot(2, 2, 2);
    hold on;
    plot(Res_final.Wing.Yst, Res_final.Wing.cdi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Wing CDi');
    plot(Res_final.Section.Y, Res_final.Section.Cd, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Section CD');
    plot(wingY_final, CD_total_final, 'g-', 'LineWidth', 2, 'DisplayName', 'Total CD');
    xlabel('Spanwise Location');
    ylabel('Drag Coefficient');
    title(sprintf('Drag Components - %s', Res_final.designName));
    legend show;
    grid on;
    hold off;
    
    % Subplot 3: Overlay of Wing CDi comparison
    subplot(2, 2, 3);
    hold on;
    plot(Res_initial.Wing.Yst, Res_initial.Wing.cdi, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Wing CDi - %s', Res_initial.designName));
    plot(Res_final.Wing.Yst, Res_final.Wing.cdi, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Wing CDi - %s', Res_final.designName));
    xlabel('Spanwise Location');
    ylabel('Wing CDi');
    title('Wing Induced Drag Comparison');
    legend show;
    grid on;
    hold off;
    
    % Subplot 4: Overlay of Total CD comparison
    subplot(2, 2, 4);
    hold on;
    plot(wingY_initial, CD_total_initial, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Total CD - %s', Res_initial.designName));
    plot(wingY_final, CD_total_final, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Total CD - %s', Res_final.designName));
    xlabel('Spanwise Location');
    ylabel('Total Drag Coefficient');
    title('Total Drag Comparison');
    legend show;
    grid on;
    hold off;
    
    % Add overall title
    sgtitle('Drag Distribution Analysis: Initial vs Final Wing Design');
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

initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;
optimizer= initializer.optimizer;
final_x_normalized = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\final.mat","x").x;
final_x = final_x_normalized.*optimizer.x0;
dvec = DesignVector().fromVector(final_x);
wingDesign = WingDesign(dvec);
optimizer.dvec = dvec;
optimizer.wingDesign = wingDesign;
optimizer.mda.wingDesign = wingDesign;

optimizer.mda.MDA_loop(Const.W_TO_max_initial,Const.W_fuel_cruise_initial,initializer.W_ZF_initial,initializer.W_AminusW_initial,initializer.V_MO_initial);

% Initial_values.optimizer.wingDesign
% optimizer.wingDesign
% After your existing code, call the comparison function
% Extract wing designs
% wingDesign_initial = initial_values.optimizer.wingDesign;
% wingDesign_final = optimizer.wingDesign;
% 
% % % Calculate results for initial design
% % Res_initial = runQ3DAnalysis(wingDesign_initial, initializer, 'Initial Design',initial_values.optimizer.mda.W_TO_max);
% % 
% % % Calculate results for final design
% % Res_final = runQ3DAnalysis(wingDesign_final, initializer, 'Final Design',optimizer.mda.W_TO_max);
% 
% plotOverlappingLift(Res_initial, Res_final)
