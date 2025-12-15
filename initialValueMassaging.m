function plotWingWeightSweep(mda, Const)

    % Sweep range for W_TO_max
    W0 = Const.W_TO_max_initial;
    sweep_range = linspace(W0 - 50000, W0 + 50000, 10);

    W_wing_vals = nan(size(sweep_range));

    for k = 1:length(sweep_range)
        W_TO = sweep_range(k);

        try
            % 1. recompute loads for this W_TO
            [lift_distribution, moment_distribution] = ...
                mda.loadsFunc(W_TO, Const.W_fuel_initial);

            % 2. recompute structures for this W_TO & new loads
            [W_TO_max_out, W_ZF_out, W_wing_out] = ...
                mda.structuresFunc(lift_distribution, ...
                                   moment_distribution, ...
                                   W_TO, ...
                                   Const.W_ZF_initial);

            % 3. store output unless invalid
            disp("Wing weight at W_TO_max "+ string(W_TO)+" is equal to "+W_wing_out);
            if ~isnan(W_wing_out)
                W_wing_vals(k) = W_wing_out;
            end

        catch
            % If loads or structure calc fails → leave NaN
            W_wing_vals(k) = NaN;
        end
    end

    % -------- PLOT --------
    figure; hold on; grid on;
    plot(sweep_range, W_wing_vals, 'LineWidth', 2);

    xlabel('W_{TO\_max} (swept)');
    ylabel('W_{wing}');
    title('Wing Weight Sensitivity to W\_{TO\_max}');
    set(gca, 'FontSize', 12);

end
function plotWingWeightSweep_ZF(mda, Const)

    % Fixed value of W_TO for all evaluations
    W_TO_fixed = Const.W_TO_max_initial;

    % Sweep range for W_ZF (±50,000 around initial)
    Z0 = Const.W_ZF_initial;
    sweep_range = linspace(Z0 - 50000, Z0 + 50000, 10);

    % storage
    W_wing_vals = nan(size(sweep_range));

    for k = 1:length(sweep_range)
        W_ZF = sweep_range(k);

        try
            % 1. loads do NOT depend on W_ZF, so use fixed W_TO
            [lift_distribution, moment_distribution] = ...
                mda.loadsFunc(W_TO_fixed, Const.W_fuel_initial);

            % 2. structures uses swept W_ZF
            [W_TO_max_out, W_ZF_out, W_wing_out] = ...
                mda.structuresFunc(lift_distribution, ...
                                   moment_distribution, ...
                                   W_TO_fixed, ...
                                   W_ZF);

            % 3. print status
            disp("Wing weight at W_ZF = " + string(W_ZF) + ...
                 " is " + string(W_wing_out));

            % 4. store wing weight (NaN will be kept for breaks)
            W_wing_vals(k) = W_wing_out;

        catch
            % leave as NaN
            W_wing_vals(k) = NaN;
        end
    end

    % -------- PLOT --------
    figure; hold on; grid on;
    plot(sweep_range, W_wing_vals, 'LineWidth', 2);

    xlabel('W_{ZF} (swept)');
    ylabel('W_{wing}');
    title('Wing Weight Sensitivity to W_{ZF}');
    set(gca, 'FontSize', 12);

end

clear all; close all; clc;

dvec = DesignVector();
wingDesign = WingDesign(dvec);
mda = MDA(wingDesign);

% plotWingWeightSweep(mda, Const);
plotWingWeightSweep_ZF(mda, Const);
