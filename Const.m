classdef Const
    properties (Constant)
        n_max = 2.5
        Mcr_ref = 0.8; 
        hcr_ref = 11673.84;    % meters
        Vcr_ref = 236.05;
        rho_fuel = 0.81715e3; % kg/m3
        f_tank = 0.93;
        W_TO_max_initial = 122470; %kg
        C_D_AnoW = 0.0248; % This was run for inviscid sim, fix update whan viscus works
        W_fuel_initial = 27220%Calculated using MTOW - W_ZF %33785; % kg
        W_ZF_initial = 95250%122470-33785; % kg
        W_wing_initial = 1.3409e+04;
        W_AminusW = 95250 - 1.3409e+04;
        LD_initial = 16;
        CD_initial = 0.0335;
        CL_initial = 0.7103;
        CT_bar = 1.8639e-4;
        V_MO_ref = 253.758;
        rho_ref = 0.327234791343548;
        drag_fus_initial = 1.6251e+04;
        q_initial = 9.1170e+03;
        internal_tank_volume =26.1383;
        S_ref = 81.8124; %[m^2]
        c_root_ref = 8.2;
        c_kink_ref = 4.66;
        c_tip_ref = 1.73;
        b_inboard_ref = 6.5;
        b_outboard_ref = 12.525
        MAC_ref = 5.04;
        LE_sweep_ref = 0.510488 %rad
        alpha_ref = -2.8384 %deg
        CD_i_ref = 0.0154725 % Induced drag
        CD_p_ref = 0.0147834 % Profile Drag
        % prop_eff, , Range
    end
end