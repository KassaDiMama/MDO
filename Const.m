classdef Const
    properties (Constant)
        n_max = 2.5
        Mcr_ref = 0.8; 
        hcr_ref = 11673.84;    % meters
        Vcr_ref = 236.05;
        rho_fuel = 0.81715e3; % kg/m3
        f_tank = 0.93;
        W_TO_max_initial = 122470; %kg
        W_fuel_initial = 33785; % kg
        W_ZF_initial = 95250%122470-33785; % kg
        W_wing_initial = 5.4582e+03;
        W_AminusW = 88685 - 5.4582e+03;
        LD_initial = 16;
        CD_initial = 0.030824958721607;
        CL_initial = 0.6853;
        CT_bar = 1.8639e-4;
    end
end