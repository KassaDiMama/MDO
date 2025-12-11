classdef Const
    properties (Constant)
        n_max = 2.5
        Mcr_ref = 0.8; 
        hcr_ref = 11673.84;    % meters
        rho_fuel = 0.81715e3; % kg/m3
        f_tank = 0.93;
        W_TO_max_initial = 122470; %kg
        W_fuel_initial = 33785; % kg
        W_ZF_initial = 122470-33785; % kg
        W_AminusW = 88685 - 5047.38;
    end
end