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
createLoadingFile(Res,"test",AC.Aero.rho,AC.Aero.V)


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
mda = MDA(wingDesign);
[lift_distribution, moment_distribution] = mda.loadsFunc(Const.W_TO_max_initial,Const.W_fuel_initial);

mda.structuresFunc(lift_distribution,moment_distribution,Const.W_TO_max_initial,Const.W_ZF_initial);