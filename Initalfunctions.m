function [rho,a, T] = isa_func(h)
% ISA_DENSITY_SIMPLE - Calculate ISA air density using piecewise linear approximation
% Input: h - altitude in meters (0 to 86,000 m)
% Output: rho - air density in kg/m³
    gamma = 1.4;
    % Convert altitude to kilometers for calculations
    h_km = h / 1000;
    
    if h_km < 11
        % Troposphere (0-11 km)
        T0 = 288.15;      % Sea level temperature (K)
        p0 = 101325;      % Sea level pressure (Pa)
        L = -6.5;         % Temperature lapse rate (K/km)
        
        T = T0 + L * h_km;
        p = p0 * (T / T0).^(-9.80665 / (L / 1000 * 287.05));
        
    elseif h_km < 20
        % Lower Stratosphere (11-20 km)
        T = 216.65;       % Constant temperature (K)
        p11 = 22632;      % Pressure at 11 km (Pa)
        
        p = p11 * exp(-9.80665 * (h_km - 11) * 1000 / (287.05 * T));
        
    
        
    else
        error('Altitude must be between 0 and 32,000 m');
    end
    
    % Calculate density using ideal gas law
    R = 287.05;          % Specific gas constant for air (J/(kg·K))
    rho = p / (R * T); a = sqrt(gamma*R*T);
    
end

[rho,a, T] = isa_func(11887.2)

function MAC = mac_func(c_r,c_t,c_kink,b_outboard)
    % Calculate the MAC of the tapered and kinked wing
    tap1 = c_kink/c_r; % taper ratio before kink
    tap2 = c_t/c_kink; %taper ratio after kink

    mac1 = 2/3 *c_r*(1+tap1+tap1^2)/(1+tap1);
    mac2 = 2/3* c_kink*(1+tap2+tap2^2)/(1+tap2);   
    S1 = 6.5/2*(c_r+c_kink);
    S2 = b_outboard/2*(c_kink+c_t);

    MAC = (S1*mac1+S2*mac2)/(S1+S2);
end

  %MAC = mac_func(8.2,1.73,28.55,12.525)

function Re = re_func(h,V,MAC)
  % Calculate Reynolds number using the formula Re = (rho * V * h) / mu
    [rho,~,T] = isa_func(h);
    Tref = 273;
    S = 110.4;
    mu0 = 1.716*10^(-5);
    mu = mu0*((T/ Tref)^1.5) * ((Tref+S)/(T+S));
    % Calculate the Reynolds number
    Re = (rho * V * MAC) / mu;
end

%Re = re_func(11887.2, 237.159, 5.0417)

function LE_sweep = le_sweep_func(c_root,c_kink)
    % Calculate leading edge sweep angle in degrees
    b = 6.5;  % Default kink location at 6.5m from root
    
    LE_sweep = atan((c_root - c_kink) / b); % Leading edge sweep angle in radians
    LE_sweep = rad2deg(LE_sweep); % Convert to degrees
end

function airfoil_x_coord = aifoil_X_coord_func(c_root,c_kink,b)
    % Calculate x coordniate of the airfoil section
    if nargin < 2 || isempty(b)
        b = 6.5;  % Default kink location at 6.5m from root
    end
    
    LE_sweep = le_sweep_func(c_root,c_kink); % Leading edge sweep angle in degrees
    % Calculate x-coordinate
    airfoil_x_coord = tan(deg2rad(LE_sweep)) * b;
end

airfoil_x_coord = aifoil_X_coord_func(8.2,4.66,19.025)

function airfoil_z_coord = aifoil_Z_coord_func(b)
    % Calculate x coordniate of the airfoil section
    
    % Calculate x-coordinate
    airfoil_z_coord = tan(deg2rad(5)) * b;
end

airfoil_z_coord = aifoil_Z_coord_func(19.025)

function S = SA_Wing_func(c_r,c_t,c_kink, b_outboard)
    % Calculate the surface area of the tapered and kinked wing
    S1 = 6.5/2*(c_r+c_kink);
    S2 = b_outboard/2*(c_kink+c_t);

    S = S1+S2;
end



function Cl = liftcoef_func(WTo_max,Wfuel,h,V,c_r,c_t,c_kink, b_outboard)
%Calculate the required lift coeficient of the aircraft at cruise
S = SA_Wing_func(c_r,c_t,c_kink, b_outboard); %sing suraface area
[rho,~,~] = isa_func(h);


L = sqrt(WTo_max*(WTo_max-Wfuel));
Cl = L/(1/2*rho*V^2*S);
end

%Cl = liftcoef_func(122470,33785, 11887.2,237.159,8.2,1.73, 28.55, 12.525)
