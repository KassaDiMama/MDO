% designVector = DesignVector();
% const = Const();
% wingDesign = WingDesign(designVector);
% 

clear all
close all
clc


%% Test create loading

% Wing planform geometry 
%                x    y     z   chord(m)    twist angle (deg) 
AC.Wing.Geom = [0     0     0     3.5         0;
                0.9  14.5   0     1.4         0];

% Wing incidence angle (degree)
AC.Wing.inc  = 0;   
            
            
% Airfoil coefficients input matrix
%                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
AC.Wing.Airfoils   = [0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797;
                      0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797];
                  
AC.Wing.eta = [0;1];  % Spanwise location of the airfoil sections

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
% catch error
%     CD = inf;
% end

%Change in line 12 the x value for exercise 3
t=toc;
% 0.2252    0.0904    0.2445    0.1314    0.4253   -0.2413   -0.0557   -0.3071  -0.2166    0.4077;
createLoadingFile(Res,"test",AC.Aero.rho,AC.Aero.V);


%% test emwet_wrapper
% 
% W_to = 52390;
% W_zf = 46720;
% 
% % emwet_wrapper(wingDesign,const, "test",W_to,W_zf)
% 
% 

%% test create airfoil
N1 = 0.5;
N2 = 1;
AU = [0.4 0.5 0.5 0.5 0.5 0.1];
% AL = [-1 -1.5 -0.9 -0.5 -0.5 -1];
AL=-AU;

[t_upper,y_upper,t_lower, y_lower,] = createAirfoilDat(N1,N2,AU,AL,"test");

% Build one continuous contour (upper forward, lower reversed)
t = [t_upper(:); flipud(t_lower(:))];
y = [y_upper(:); flipud(y_lower(:))];

% Optional: close the loop by repeating the first point
t(end+1) = t(1);
y(end+1) = y(1);

plot(t, y, '-k');
axis equal;
xlabel('t'); ylabel('y');
grid on;

%% Testing volume wing
clear all

dvec = DesignVector();
wingDesign = WingDesign(dvec);

%% Calculating Initial W_AminusW
clear all
close all
clc
dvec = DesignVector();
wingDesign = WingDesign(dvec);
mda = MDA(wingDesign,Const.W_TO_max_initial,Const.W_ZF_initial);
[lift_distribution, moment_distribution] = mda.loadsFunc(Const.W_TO_max_initial);

[W_TO_max_out, W_ZF_out, W_wing_out] = mda.structuresFunc(lift_distribution, moment_distribution, Const.W_TO_max_initial, Const.W_ZF_initial);
%% MDA LOOP TEST
clear all
close all
clc
dvec = DesignVector();
wingDesign = WingDesign(dvec);
mda = MDA(wingDesign);
mda.MDA_loop(Const.W_TO_max_initial,Const.W_fuel_initial,wingDesign.W_fuel);
%% Check Loading
clear all
close all
clc

dvec = DesignVector();
wingDesign = WingDesign(dvec);
mda = MDA(wingDesign);
[lift_distribution, moment_distribution] = mda.loadsFunc(Const.W_TO_max_initial,wingDesign.W_fuel);


% Plot lift and moment distributions
figure

subplot(2,1,1)
plot(lift_distribution.y, lift_distribution.L, 'LineWidth', 1.5)
grid on
xlabel('Spanwise Location y [m]')
ylabel('Lift')
title('Lift Distribution')

subplot(2,1,2)
plot(moment_distribution.y, moment_distribution.M, 'LineWidth', 1.5)
grid on
xlabel('Spanwise Location y [m]')
ylabel('Moment')
title('Moment Distribution')
%% Plot wing
dvec = DesignVector();
% dvec.LE_sweep = 20/180*pi;
wingDesign = WingDesign(dvec);

x_root = wingDesign.x_root;
x_kink = wingDesign.x_kink;
x_tip  = wingDesign.x_tip;

y_root = wingDesign.y_root;
y_kink = wingDesign.y_kink;
y_tip  = wingDesign.y_tip;

c_root = wingDesign.c_root;
c_kink = wingDesign.c_kink;
c_tip  = wingDesign.c_tip;

% Calculate trailing edge coordinates (top view)
x_te_root = x_root + c_root;
x_te_kink = x_kink + c_kink;
x_te_tip  = x_tip  + c_tip;

% Leading edge coordinates (top view)
LE_x = [x_root, x_kink, x_tip];
LE_y = [y_root, y_kink, y_tip];

% Trailing edge coordinates (top view)
TE_x = [x_te_root, x_te_kink, x_te_tip];
TE_y = [y_root, y_kink, y_tip];

% Combine for plotting the wing outline
wing_x = [LE_x, fliplr(TE_x)];
wing_y = [LE_y, fliplr(TE_y)];

% Plot
figure;
fill(wing_x, wing_y, [0.6 0.8 1]); % Wing shape filled with color
hold on;
plot(LE_x, LE_y, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor','k'); % Leading edge
plot(TE_x, TE_y, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor','r'); % Trailing edge
axis equal
xlabel('x (m)');
ylabel('y (m)');
title('Top-down view of wing planform');
grid on;
legend('Wing','Leading Edge','Trailing Edge');
% %% Fuselage drag
% clear all
% 
% designVector = DesignVector();
% const = Const();
% wingDesign = WingDesign(designVector);
% 
% function obj = MDA(wingDesign)
%             arguments
%                 wingDesign WingDesign
%             end
%             wingDesign = wingDesign ;
% end
% 
% obj = MDA(wingDesign)
% % Wing planform geometry 
% %               x    y     z   chord(m)    twist angle (deg) 
% AC.Wing.Geom = [obj.wingDesign.x_root     obj.wingDesign.y_root     obj.wingDesign.z_root     obj.wingDesign.c_root         obj.wingDesign.twist(1)
%                 obj.wingDesign.x_kink     obj.wingDesign.y_kink     obj.wingDesign.z_kink     obj.wingDesign.c_kink         obj.wingDesign.twist(2)
%                 obj.wingDesign.x_tip     obj.wingDesign.y_tip     obj.wingDesign.z_tip     obj.wingDesign.c_tip        obj.wingDesign.twist(3)];
% % AC.Wing.Geom = [0     0     0     3.5         0;
% %     0.9  14.5   0     1.4         0
% %     2*0.9  2*14.5   0     1.4         0];
% % Wing incidence angle (degree)
% AC.Wing.inc  = obj.wingDesign.incidence;   
% 
% 
% % Airfoil coefficients input matrix
% %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
% AC.Wing.Airfoils   = [obj.wingDesign.AU obj.wingDesign.AL;
%                       obj.wingDesign.AU obj.wingDesign.AL];
% 
% %AC.Wing.eta = [obj.wingDesign.y_root/obj.wingDesign.b_half;obj.wingDesign.y_kink/obj.wingDesign.b_half];  % Spanwise location of the airfoil sections
% AC.Wing.eta = [0;1];
% % Viscous vs inviscid
% AC.Visc  = 0;              % 0 for inviscid and 1 for viscous analysis
% AC.Aero.MaxIterIndex = 150;
% % Flight Condition
% AC.Aero.V     = obj.wingDesign.V;          % flight speed (m/s)
% AC.Aero.rho   = obj.wingDesign.rho;         % air density  (kg/m3)
% AC.Aero.alt   = obj.wingDesign.hcr;             % flight altitude (m)
% AC.Aero.Re    = obj.wingDesign.Re;        % reynolds number (bqased on mean aerodynamic chord)
% AC.Aero.M     = obj.wingDesign.Mcr;          % flight Mach number 
% AC.Aero.CL    = obj.wingDesign.liftcoef_func(Const.W_TO_max_initial,const.W_fuel_initial);          % lift coefficient - comment this line to run the code for given alpha%
% % AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 
% 
% % tic
% 
% % try 
% % disp("Starting Q3D");
% Res = Q3D_solver(AC);
% 
% AC.Aero.M
% AC.Aero.V
% Cdw = drag_estimation(Res,AC.Visc,true)
% [CDi_total, CDv_total] = Drag_coeff_from_spanwise(Res, wingDesign.S);
% fprintf("induced drag equals: %g\n",CDi_total);
% fprintf("profile drag equals: %g\n",CDv_total);
% % Cd_AnoW = Res.CLwing/16 - Cdw


%% Calculate initial lift and drag
clear all
close all
clc

dvec = DesignVector();
optimizer = Optimizer(dvec);
[CL_wing, CD_wing] = optimizer.calcCL_CD(Const.W_TO_max_initial,optimizer.wingDesign.W_fuel);

%% Calculate range
clear all
close all
clc

dvec = DesignVector();
initializer = Initializer(dvec);
optimizer = initializer.optimizer;
x = dvec.toVector();
range = optimizer.objective_wrapper(x./x); % in meters
fprintf("Initial range equals: %g km\n",-range/1000);


%% Reference aircraft values
clear all

dvec = DesignVector();
wingDesign = WingDesign(dvec);
const = Const();

fprintf("Wing surface area S equals: %g\n",wingDesign.S);
fprintf("Wing MAC equals: %g\n",wingDesign.MAC);
fprintf("Wing LE Sweep equals: %g\n",wingDesign.LE_sweep);

%% Run the optimization
clear all
close all
clc

echo all off
dvec = DesignVector();
initializer = Initializer(dvec);
save("initializer"+datestr(now,'yyyy-mm-dd_HH-MM-SS')+".mat", 'initializer');
optimizer = initializer.optimizer;
msg = [
    "---------------------------------"
    "---------------------------------"
    "Starting new run at " + string(datestr(now, 'yyyy-mm-dd HH:MM:SS'))
];
fname = "run_" + datestr(now,'yyyy-mm-dd_HH-MM-SS') + ".txt";
diary(fname)
% diary on
% echo Optimizer.m on
% echo MDA.m on
logMessage(msg, "log.file")
optimizer.start();
%% Test Initializer
clear all
close all
clc
dvec = DesignVector();
initializer = Initializer(dvec);

%% See result
clear all
close all
clc

function plotWing(ax, wingDesign, varargin)
    % plotWing  Plot wing planform into specified axes
    %
    % Inputs:
    %   ax         - axes handle
    %   wingDesign - struct or object with wing geometry fields
    %
    % Name-value pairs (optional):
    %   'Color'    - RGB triplet or color char (default: [0.6 0.8 1])
    %   'Label'    - Legend label (default: 'Wing')
    %   'Alpha'    - Face transparency (default: 0.6)

    % --- Defaults ---
    p = inputParser;
    addParameter(p,'Color',[0.6 0.8 1]);
    addParameter(p,'Label','Wing');
    addParameter(p,'Alpha',0.6);
    parse(p,varargin{:});

    wingColor = p.Results.Color;
    wingLabel = p.Results.Label;
    alphaVal  = p.Results.Alpha;

    % --- Extract geometry ---
    x_root = wingDesign.x_root;
    x_kink = wingDesign.x_kink;
    x_tip  = wingDesign.x_tip;
    
    y_root = wingDesign.y_root;
    y_kink = wingDesign.y_kink;
    y_tip  = wingDesign.y_tip;
    
    c_root = wingDesign.c_root;
    c_kink = wingDesign.c_kink;
    c_tip  = wingDesign.c_tip;
    
    % --- Trailing edge ---
    x_te_root = x_root + c_root;
    x_te_kink = x_kink + c_kink;
    x_te_tip  = x_tip  + c_tip;
    
    % --- Leading & trailing edges ---
    LE_x = [x_root, x_kink, x_tip];
    LE_y = [y_root, y_kink, y_tip];
    
    TE_x = [x_te_root, x_te_kink, x_te_tip];
    TE_y = [y_root, y_kink, y_tip];
    
    % --- Wing outline ---
    wing_x = [LE_x, fliplr(TE_x)];
    wing_y = [LE_y, fliplr(TE_y)];
    
    % --- Plot ---
    hWing = fill(ax, wing_x, wing_y, wingColor, ...
        'FaceAlpha', alphaVal, ...
        'EdgeColor', 'none', ...
        'DisplayName', wingLabel);
    hold(ax, 'on');

    plot(ax, LE_x, LE_y, 'k-', 'LineWidth', 1.5, ...
        'HandleVisibility','off');

    plot(ax, TE_x, TE_y, 'k--', 'LineWidth', 1.5, ...
        'HandleVisibility','off');

    axis(ax, 'equal');
    grid(ax, 'on');

    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    title(ax, 'Top-down view of wing planform');

    legend(ax, 'show');
end


initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;

fminconresults = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\2025-12-24_20-15-52fmincon_results.mat");

x_opt_normalized = fminconresults.x_opt;
x_opt = x_opt_normalized .* initializer.optimizer.x0;

dvec = DesignVector();
dvec = dvec.fromVector(x_opt);
% dvec.LE_sweep = 20/180*pi;
wingDesign_new = WingDesign(dvec);


figure
ax = axes;
hold(ax,'on')



plotWing(ax, wingDesign_new, ...
    'Color',[1 0.4 0.4], ...
    'Label','Optimized Wing');
fail_wingDesign = load("emwet_fail.mat").ans;
% plotWing(ax, fail_wingDesign, ...
%     'Color',[1 0.4 0.4], ...
%     'Label','Optimized Wing');
plotWing(ax, initializer.optimizer.wingDesign, ...
    'Color',[0.6 0.8 1], ...
    'Label','Baseline Wing');

%% Test taper ratio c_kink
clear all
close all
clc
dvec = DesignVector();
initializer = Initializer(dvec);
TR_new = initializer.getTRbounds();

dvec.TR= TR_new;
dvec.AR = Const.AR_upper_bound;
dvec.LE_sweep = Const.LE_sweep_upper_bound/180*pi;
dvec.b_half = Const.b_half_upper_bound;
% dvec.AR= dvec.AR*1.1;
wingDesign = WingDesign(dvec);
% fprintf("Changing AR value in design vector.\n");
% fprintf("c_kink equals: %g\n", wingDesign.c_kink);
% fprintf("c_tip equals: %g\n", wingDesign.c_tip);
% ratio = wingDesign.c_kink / wingDesign.c_tip;
% fprintf("Ratio of c_kink to c_tip equals: %g\n", ratio);
% 
% 
% dvec.TR= TR_new;
% dvec.AR= dvec.AR/1.1;
% dvec.LE_sweep= dvec.LE_sweep*1.1;
% wingDesign = WingDesign(dvec);
% fprintf("Changing LE_sweep value in design vector.\n");
% fprintf("c_kink equals: %g\n", wingDesign.c_kink);
% fprintf("c_tip equals: %g\n", wingDesign.c_tip);
% ratio = wingDesign.c_kink / wingDesign.c_tip;
% fprintf("Ratio of c_kink to c_tip equals: %g\n", ratio);
% 
% dvec.TR= TR_new;
% dvec.LE_sweep= dvec.LE_sweep/1.1;
% dvec.b_half= dvec.b_half*0.9;
% wingDesign = WingDesign(dvec);
% fprintf("Changing b_half value in design vector.\n");
% fprintf("c_kink equals: %g\n", wingDesign.c_kink);
% fprintf("c_tip equals: %g\n", wingDesign.c_tip);
% ratio = wingDesign.c_kink / wingDesign.c_tip;
% fprintf("Ratio of c_kink to c_tip equals: %g\n", ratio);
% 
% dvec.TR= TR_new*1.1;
% dvec.b_half= dvec.b_half/0.9;
% wingDesign = WingDesign(dvec);
% fprintf("Changing TR value in design vector.\n");
% fprintf("c_kink equals: %g\n", wingDesign.c_kink);
% fprintf("c_tip equals: %g\n", wingDesign.c_tip);
% ratio = wingDesign.c_kink / wingDesign.c_tip;
% fprintf("Ratio of c_kink to c_tip equals: %g\n", ratio);


x_root = wingDesign.x_root;
x_kink = wingDesign.x_kink;
x_tip  = wingDesign.x_tip;

y_root = wingDesign.y_root;
y_kink = wingDesign.y_kink;
y_tip  = wingDesign.y_tip;

c_root = wingDesign.c_root;
c_kink = wingDesign.c_kink;
c_tip  = wingDesign.c_tip;

% Calculate trailing edge coordinates (top view)
x_te_root = x_root + c_root;
x_te_kink = x_kink + c_kink;
x_te_tip  = x_tip  + c_tip;

% Leading edge coordinates (top view)
LE_x = [x_root, x_kink, x_tip];
LE_y = [y_root, y_kink, y_tip];

% Trailing edge coordinates (top view)
TE_x = [x_te_root, x_te_kink, x_te_tip];
TE_y = [y_root, y_kink, y_tip];

% Combine for plotting the wing outline
wing_x = [LE_x, fliplr(TE_x)];
wing_y = [LE_y, fliplr(TE_y)];

% Plot
figure;
fill(wing_x, wing_y, [0.6 0.8 1]); % Wing shape filled with color
hold on;
plot(LE_x, LE_y, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor','k'); % Leading edge
plot(TE_x, TE_y, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor','r'); % Trailing edge
axis equal
xlabel('x (m)');
ylabel('y (m)');
title('Top-down view of wing planform');
grid on;
legend('Wing','Leading Edge','Trailing Edge');

%% Drag plots
clear all
initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;
optimizer= initializer.optimizer;
final_x_normalized = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\final.mat","x").x;
final_x = final_x_normalized.*optimizer.x0;
dvec = DesignVector().fromVector(final_x);
wingDesign = WingDesign(dvec);
optimizer.dvec = dvec;
optimizer.wingDesign = wingDesign;
optimizer.mda.wingDesign = wingDesign;
const = Const();

optimizer.mda.MDA_loop(Const.W_TO_max_initial,Const.W_fuel_cruise_initial,initializer.W_ZF_initial,initializer.W_AminusW_initial,initializer.V_MO_initial);
% 
% function Result = Q3D(wingDesign)
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
                  
%AC.Wing.eta = [wingDesign.y_root/wingDesign.b_half;wingDesign.y_kink/wingDesign.b_half;wingDesign.y_tip/wingDesign.b_half];  % Spanwise location of the airfoil sections
AC.Wing.eta = [0;1];
% Viscous vs inviscid
AC.Visc  = 0;              % 0 for inviscid and 1 for viscous analysis
AC.Aero.MaxIterIndex = 150;
% Flight Condition
[rho_dont_use,a,T_dont_use] = wingDesign.isa_func();
Mcritical = initializer.V_MO_initial/a;
Re_corrected = wingDesign.Re/initializer.V_MO_initial*wingDesign.V;
AC.Aero.V     = initializer.V_MO_initial;            % flight speed (m/s)
AC.Aero.rho   = wingDesign.rho;         % air density  (kg/m3)
AC.Aero.alt   = wingDesign.hcr;             % flight altitude (m)
AC.Aero.Re    = Re_corrected;        % reynolds number (bqased on mean aerodynamic chord)
AC.Aero.M     = Mcritical;           % flight Mach number 
AC.Aero.CL    = wingDesign.calculateCL_cruise(const.W_TO_max_initial,initializer.V_MO_initial);          % lift coefficient - comment this line to run the code for given alpha%
% AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 

% tic

% try 
% disp("Starting Q3D");
Res = Q3D_solver(AC)



%% Example
clear all
close all
clc
dvec1 = DesignVector();
initial_values = Initializer(dvec1);

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

%% Example

final_x = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\final.mat","x").x;
dvec = DesignVector().fromVector(final_x);
wingDesign = WingDesign(dvec);
initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;
optimizer= Optimizer(dvec,wingDesign,initializer);
optimizer.mda.MDA_loop();
optimizer.mda.loadsFunc(optimizer.mda.W_TO_max,initializer.V_MO_initial)

%% Example
clear all
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
W_fuel = optimizer.mda.W_TO_max-optimizer.mda.W_ZF;
[CL_wing, CD_wing]=optimizer.calcCL_CD(optimizer.mda.W_TO_max,W_fuel);

%% Calculating Initial Values For Report

initializer = load("fmincon_2026-01-07_11-55-38\initializer2026-01-07_11-51-48.mat").initializer;
W_wing = initializer.W_wing_initial;
emission = Const.W_fuel_cruise_initial*3.16;
wing_tank_volume = initializer.wing_tank_volume_initial;
eta = initializer.optimizer.performanceFunction();
CL_wing = initializer.CL_initial;
[alpha_cruise, CD_wing, CD_wing_induced] = initializer.optimizer.calculateAoACruise(Const.W_TO_max_initial, Const.W_fuel_cruise_initial);
wingless_drag_force = initializer.drag_fus_initial;
cd_fus_init=wingless_drag_force/initializer.optimizer.calculateDesignDynamicPressure()/initializer.S_initial;
W_a_min_w = initializer.W_AminusW_initial;
fprintf('W_wing: %f\n', W_wing);
fprintf('Emission: %f\n', emission);
fprintf('Wing Tank Volume: %f\n', wing_tank_volume);
fprintf('Efficiency (eta): %f\n', eta);
fprintf('CL_wing: %f\n', CL_wing);
fprintf('Alpha Cruise: %f\n', alpha_cruise);
fprintf('CD_wing: %f\n', CD_wing);
fprintf('CD_wing Induced: %f\n', CD_wing_induced);
fprintf('Wingless Drag Force: %f\n', wingless_drag_force);
fprintf('CD_fus_init: %f\n', cd_fus_init);
fprintf('W_a_min_w: %f\n', W_a_min_w);
CT = Const.CT_bar/eta;
fprintf('CT: %f\n', CT);

%% Calculate Optimized
clear all
initializer = load("fmincon_2026-01-07_11-55-38\initializer2026-01-07_11-51-48.mat").initializer;
optimizer= initializer.optimizer;
final_x_normalized = load("fmincon_2026-01-07_11-55-38\final.mat","x").x;
final_x = final_x_normalized.*optimizer.x0;
dvec = DesignVector().fromVector(final_x);
wingDesign = WingDesign(dvec);
optimizer.dvec = dvec;
optimizer.wingDesign = wingDesign;
optimizer.mda.wingDesign = wingDesign;

optimizer.mda.MDA_loop(Const.W_TO_max_initial,Const.W_fuel_cruise_initial,initializer.W_ZF_initial,initializer.W_AminusW_initial,initializer.V_MO_initial);
W_fuel = optimizer.mda.W_TO_max-optimizer.mda.W_ZF;
[CL_wing, CD_wing]=optimizer.calcCL_CD(optimizer.mda.W_TO_max,W_fuel);
W_TO = optimizer.mda.W_TO_max;
[lift_dist, moment_dist] = optimizer.mda.loadsFunc(W_TO,initializer.V_MO_initial);
W_wing = optimizer.mda.structuresFunc(lift_dist,moment_dist,W_TO,optimizer.mda.W_ZF);
W_co2 = W_fuel * 3.16;
FuelVolume = W_fuel / 0.81715e3;
tank_volume = optimizer.wingDesign.calculateWingTankVolume();
dynamic_pressure=optimizer.calculateDesignDynamicPressure();
eta=optimizer.performanceFunction();
CL_cr = optimizer.wingDesign.calculateCL_cruise(W_TO,W_fuel);
V = optimizer.wingDesign.V;
hcr = optimizer.wingDesign.hcr;
mcr = optimizer.wingDesign.Mcr;
Re_cr = optimizer.wingDesign.Re;
[alpha_cruise, CD_wing, CD_wing_induced] = optimizer.calculateAoACruise(Const.W_TO_max_initial, Const.W_fuel_cruise_initial);
L_over_D = optimizer.aerodynamicsFunc(W_TO,W_fuel);
W_a_min_w = W_TO - W_fuel-W_wing;
S = optimizer.wingDesign.S;
MAC = optimizer.wingDesign.MAC;
wing_loading = W_TO/S;
AR = optimizer.wingDesign.AR;
LE_sweep = optimizer.wingDesign.LE_sweep;
b_inboard = optimizer.wingDesign.b_inboard;
b_outboard = optimizer.wingDesign.b_outboard;
c_root = optimizer.wingDesign.c_root;
c_kink = optimizer.wingDesign.c_kink;
c_tip = optimizer.wingDesign.c_tip;

range = optimizer.objectiveFunc(optimizer.wingDesign.W_fuel, optimizer.mda.W_TO_max, L_over_D, eta);
objective = -(range/optimizer.initializer.range_initial);

wing_loading_constraint_value = (optimizer.mda.W_TO_max/optimizer.wingDesign.S  - Const.W_TO_max_initial/optimizer.initializer.S_initial)/(Const.W_TO_max_initial/optimizer.initializer.S_initial); % Wing loading constraint  


fprintf('W_TO: %f\n', W_TO);
fprintf('W_fuel: %f\n', W_fuel);
fprintf('W_wing: %f\n', W_wing);
fprintf('W_co2: %f\n', W_co2);
fprintf('Fuel Volume: %f\n', FuelVolume);
fprintf('Tank Volume: %f\n', tank_volume);
fprintf('Dynamic Pressure: %f\n', dynamic_pressure);
fprintf('Efficiency (eta): %f\n', eta);
fprintf('CL_cr: %f\n', CL_cr);
fprintf('V: %f\n', V);
fprintf('hcr: %f\n', hcr);
fprintf('mcr: %f\n', mcr);
fprintf('Re_cr: %f\n', Re_cr);
fprintf('Alpha Cruise: %f\n', alpha_cruise);
fprintf('CD_wing: %f\n', CD_wing);
fprintf('CD_wing Induced: %f\n', CD_wing_induced);
fprintf('L/D Ratio: %f\n', L_over_D);
fprintf('Wing Loading: %f\n', wing_loading);
fprintf('Aspect Ratio (AR): %f\n', AR);
fprintf('Leading Edge Sweep: %f\n', LE_sweep);
fprintf('Inboard Span: %f\n', b_inboard);
fprintf('Outboard Span: %f\n', b_outboard);
fprintf('Root Chord: %f\n', c_root);
fprintf('Kink Chord: %f\n', c_kink);
fprintf('Tip Chord: %f\n', c_tip);
fprintf('Objective: %f\n', objective);
fprintf('Range: %f\n', range);
fprintf('Wing Loading Constraint Value: %f\n', wing_loading_constraint_value);

%% Analyze Convergence History
initializer = load("fmincon_2026-01-07_11-55-38\initializer2026-01-07_11-51-48.mat").initializer;
optimizer = initializer.optimizer;
iteration_history = {};
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
for iteration = 0:5
    filename = sprintf("fmincon_2026-01-07_11-55-38/iter_%05d.mat", iteration);
    data = load(filename);
    x_normalized = data.x;
    x_scaled = x_normalized .* optimizer.x0;
    dvec = DesignVector().fromVector(x_scaled);
    wingDesign = WingDesign(dvec);

    optimizer.dvec = dvec;
    optimizer.wingDesign = wingDesign;
    optimizer.mda.wingDesign = wingDesign;
    
    optimizer.mda.MDA_loop(Const.W_TO_max_initial,Const.W_fuel_cruise_initial,initializer.W_ZF_initial,initializer.W_AminusW_initial,initializer.V_MO_initial);
    eta=optimizer.performanceFunction();
    W_fuel = optimizer.mda.W_TO_max-optimizer.mda.W_ZF;
    [CL_wing, CD_wing]=optimizer.calcCL_CD(optimizer.mda.W_TO_max,W_fuel);
    W_TO = optimizer.mda.W_TO_max;
    L_over_D = optimizer.aerodynamicsFunc(W_TO,W_fuel);
    range = optimizer.objectiveFunc(optimizer.wingDesign.W_fuel, optimizer.mda.W_TO_max, L_over_D, eta);
    objective = -(range/optimizer.initializer.range_initial);
    
    wing_loading_constraint_value = (optimizer.mda.W_TO_max/optimizer.wingDesign.S  - Const.W_TO_max_initial/optimizer.initializer.S_initial)/(Const.W_TO_max_initial/optimizer.initializer.S_initial); % Wing loading constraint  
    
    W_TO = optimizer.mda.W_TO_max;
    S = optimizer.wingDesign.S;
    wing_loading = W_TO/S;

    iteration_block.iteration = iteration;
    iteration_block.range = range;
    iteration_block.objective = objective;
    iteration_block.wing_loading_constraint_value = wing_loading_constraint_value;
    iteration_block.wing_loading = wing_loading;

    N1 = 0.5;
    N2 = 1;
    CST_order = length(optimizer.wingDesign.AU) - 1;
    
    
    
    ts = linspace(0, 1, 10000);
    yu = CSTcurve(ts, optimizer.wingDesign.AU, N1, N2, CST_order);
    yl = CSTcurve(ts, optimizer.wingDesign.AL, N1, N2, CST_order);
    
    mask_u = yu(2:end-1);
    mask_l = yl(2:end-1);
    res = mask_u < mask_l;
    upperLowerOverlapFraction = sum(res)/length(mask_u);
    iteration_block.iteration = iteration;
    iteration_block.range = range;
    iteration_block.objective = objective;
    iteration_block.wing_loading_constraint_value = wing_loading_constraint_value;
    iteration_block.wing_loading = wing_loading;
    iteration_block.upperLowerOverlapFraction = upperLowerOverlapFraction;
    iteration_block.aboveUpper = sum(res);
    iteration_history{iteration+1} = iteration_block;

end

%% Plot Convergence History
% Load file
data = load("fmincon_2026-01-07_11-55-38/iteration_history2.mat");

iteration_history = data.iteration_history;

n = numel(iteration_history);

% Preallocate
iterations    = zeros(1, n);
ranges        = zeros(1, n);
wing_loadings = zeros(1, n);
objectives = zeros(1, n);
wing_loading_constraints = zeros(1, n);
aboveUpper = zeros(1,n);
upperLowerOverlapFractions = zeros(1, n); % Preallocate for upper-lower overlap fractions

for i = 1:n
    iterations(i)               = iteration_history{i}.iteration;
    ranges(i)                   = iteration_history{i}.range;
    wing_loadings(i)            = iteration_history{i}.wing_loading;
    objectives(i)               = iteration_history{i}.objective;
    wing_loading_constraints(i) = iteration_history{i}.wing_loading_constraint_value;
    aboveUpper(i)               = iteration_history{i}.aboveUpper;
    upperLowerOverlapFractions(i) = iteration_history{i}.upperLowerOverlapFraction;
end


figure;
plot(iterations, objectives, 'b-', 'LineWidth', 2)
xlabel('Iteration')
ylabel('Normalized Objective')
title('Objective Function History')
grid on

figure;
plot(iterations, wing_loading_constraints, 'b-', 'LineWidth', 2); hold on
plot(iterations, upperLowerOverlapFractions, 'r-', 'LineWidth', 2);
yline(0, '--k', 'Constraint Boundary', 'LineWidth', 1.5);
ymin = min(wing_loading_constraints)-0.01;
ymax = 0.01;
ylim([ymin ymax]);
hold off

xlabel('Iteration')
ylabel('Normalized Value')
title('Wing Loading Constraint & Upper–Lower Overlap History')
legend('Constraint','Overlap Fraction','Boundary','Location','best')
grid on

%%
display(searchdir.message)
