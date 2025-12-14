clear all
close all
clc

N1 = 0.5;
N2 = 1;
CST_order = 5;
total_points = 205;
num_per_side = (total_points+1)/2;
% AU = [0.4 0.5 0.5 0.5 0.5 0.1];
AU = 0.5 * ones(1, CST_order+1);
% AL = [-0.04 -0.04 -0.04 -0.04 -0.5 -1];
AL = 0.5 * ones(1, CST_order+1);
% AL=-AU;


try
    % airfoilData = load('withcomb135.dat');
    airfoilData = load('sc207210.dat');
    airfoil_X = airfoilData(:, 1);
    airfoil_Y = airfoilData(:, 2);
catch ME
    disp(['Error loading airfoil file: ', ME.message]);
    return; % Exit the script if the file cannot be loaded
end
ts= airfoilData(1:num_per_side, 1);
[t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1,N2,AU,AL,"test",flip(ts),CST_order);

% Build one continuous contour (upper forward, lower reversed)
t = [t_upper; t_lower];
y = [y_upper; y_lower];

% Optional: close the loop by repeating the first point
% t(end+1) = t(1);
% y(end+1) = y(1);

% plot(t, y, '-k');
% axis equal;
% xlabel('t'); ylabel('y');
% grid on;


% STEP 2: Load the Airfoil Data
% Replace 'naca0012.dat' with your actual file name
try
    airfoilData = load('sc207210.dat');
    airfoil_X = airfoilData(:, 1);
    airfoil_Y = airfoilData(:, 2);
catch ME
    disp(['Error loading airfoil file: ', ME.message]);
    return; % Exit the script if the file cannot be loaded
end

% STEP 3: Create the Plot and Add the Airfoil
figure;
plot(t, y, '-b', 'LineWidth', 1.5, 'DisplayName', 'Original Data');

% Hold the plot to add the airfoil
hold on; 

% Plot the airfoil
plot(airfoil_X, airfoil_Y, 'r-', 'LineWidth', 2, 'DisplayName', 'Airfoil Profile'); 

% Apply axis properties
axis equal; % Essential for an airfoil plot to look correct
xlabel('X (Airfoil)'); 
ylabel('Y (Airfoil)');
title('Plot with Airfoil Data');
grid on;
legend('show', 'Location', 'southwest');

hold off;

function objective = objectiveFunction(x, N1, N2, airfoilData, airfoil_Y,CST_order,num_per_side)
    AU = x(1:CST_order+1);
    AL = x(CST_order+2:2*CST_order+2);
    ts= airfoilData(1:num_per_side, 1);
    [t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1,N2,AU,AL,"test",flip(ts),CST_order);
    
    % Build one continuous contour (upper forward, lower reversed)
    y = [y_upper; y_lower];

    objective = sum((y-airfoil_Y).^2);
end 

x0 = [AU;AL];
lb = -0.1 * ones(1, CST_order+1);;
ub = -50*lb;

A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

objective = @(x) objectiveFunction(x, N1, N2, airfoilData, airfoil_Y,CST_order,num_per_side);


% Options (recommended)
options = optimoptions('fmincon',...
    'Algorithm','sqp', ...
    'Display','iter-detailed', ...
    'MaxIterations',500, ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxFunctionEvaluations', 1e6, ...
    'OptimalityTolerance',1e-10);

% Run optimization
[x_opt, fval, exitflag, output] = fmincon(objective, ...
                                          x0, A, b, Aeq, beq, lb, ub, ...
                                          nonlcon, options);

AU = x_opt(1:CST_order+1);
AL = x_opt(CST_order+2:2*CST_order+2);

ts_original= airfoilData(1:num_per_side, 1);
ts = linspace(0,1,num_per_side);
[t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1,N2,AU,AL,"test",flip(ts));

% Build one continuous contour (upper forward, lower reversed)
t = [t_upper t_lower];
y = [y_upper y_lower];

% Plot

% STEP 3: Create the Plot and Add the Airfoil
figure;
plot(t, y, '-b', 'LineWidth', 1.5, 'DisplayName', 'CST Data');

% Hold the plot to add the airfoil
hold on; 

% Plot the airfoil
plot(airfoil_X, airfoil_Y, 'r-', 'LineWidth', 2, 'DisplayName', 'Airfoil Original Profile'); 

% Apply axis properties
axis equal; % Essential for an airfoil plot to look correct
xlabel('X (Airfoil)'); 
ylabel('Y (Airfoil)');
title('solution');
grid on;
legend('show', 'Location', 'southwest');

hold off;