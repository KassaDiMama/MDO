classdef Const
    properties (Constant)
        n_max = 2.5;
        Mcr_ref = 0.8; 
        hcr_ref = 11673.84;    % meters
        % Vcr_ref = 236.05;
        rho_fuel = 0.81715e3; % kg/m3
        f_tank = 0.93;
        W_TO_max_initial = 122470; %kg
        W_fuel_cruise_initial = 33785 ;%Got from excel %33785; % kg
        LD_initial = 14;
        CT_bar = 1.8639e-4;
        internal_tank_volume =0;
        b_inboard_ref = 6.5;
        b_outboard_ref = 12.525
        LE_sweep_ref = 0.49836829593%THIS IS 25 deg qc to LE converted   + 4 / 180 *pi%0.510488 %rad THIS
        AR_ref = 7.82;
        TR_ref = 0.243;
        fuel_weight_max_ref = 34890
        airfoil_ref = 'withcomb135.dat'
        airfoil_thickness_multiplier = 1
        b_half_upper_bound = 26
        b_half_lower_bound = 15

        LE_sweep_upper_bound = 40 % deg
        LE_sweep_lower_bound = 0 % deg

        AR_upper_bound = 12
        AR_lower_bound = 6

        TR_upper_bound = 0.85
        TR_lower_bound = 0.2
        % C_D_AnoW = 0.0248; % This was run for inviscid sim, fix update whan viscus works
        
        % W_ZF_initial = 122470-33785; % kg
        % W_wing_initial = 15429.2;
        % W_AminusW = 95250 - 15429.2;
        
        % CD_initial = 0.0335;
        % CL_initial = 0.7103;
        
        % V_MO_ref = 253.758;
        % rho_ref = 0.327234791343548;
        % drag_fus_initial = 1.6251e+04;
        % q_initial = 9.1170e+03;
        
        % S_ref = 81.8124; %[m^2]
        % c_root_ref = 8.2;
        % c_kink_ref = 4.66;
        % c_tip_ref = 1.73;
        
        % MAC_ref = 5.04;
        
        % alpha_ref = -2.8384 %deg
        % CD_i_ref = 0.0154725 % Induced drag
        % CD_p_ref = 0.0147834 % Profile Drag
        
        % prop_eff, , Range
    end
end