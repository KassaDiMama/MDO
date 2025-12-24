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
%             obj.wingDesign = wingDesign ;
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
dvec.b_half=dvec.b_half+4;
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


initializer = load("initializer.mat").initializer;
fminconresults = load("fmincon_2025-12-24_18-11-55/fmincon_results.mat");

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