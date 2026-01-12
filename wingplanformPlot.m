function plotWingPlanformComparison(wingDesign1, wingDesign2, varargin)
    % plotWingPlanformComparison - Plots top-down view of two wing planforms overlapped
    %
    % Inputs:
    %   wingDesign1 - struct containing first wing geometric parameters
    %   wingDesign2 - struct containing second wing geometric parameters
    %
    % Optional Name-Value pairs:
    %   'Color1' - First wing fill color [R G B] (default: [0.6 0.8 1])
    %   'Color2' - Second wing fill color [R G B] (default: [1 0.6 0.6])
    %   'Alpha' - Face transparency (default: 0.6)
    %   'LineWidth' - Line width for edges (default: 1.5)
    %   'Label1' - Label for first wing (default: 'Wing 1')
    %   'Label2' - Label for second wing (default: 'Wing 2')
    %   'FigureHandle' - Handle to existing figure (default: creates new)
    %   'AxesHandle' - Handle to existing axes (default: creates new)
    
    % --- Defaults ---
    p = inputParser;
    addParameter(p, 'Color1', [0.6 0.8 1]);
    addParameter(p, 'Color2', [1 0.6 0.6]);
    addParameter(p, 'Alpha', 0.6);
    addParameter(p, 'LineWidth', 1.5);
    addParameter(p, 'Label1', 'Wing 1');
    addParameter(p, 'Label2', 'Wing 2');
    addParameter(p, 'FigureHandle', []);
    addParameter(p, 'AxesHandle', []);
    parse(p, varargin{:});
    
    wingColor1 = p.Results.Color1;
    wingColor2 = p.Results.Color2;
    alphaVal = p.Results.Alpha;
    lineWidth = p.Results.LineWidth;
    label1 = p.Results.Label1;
    label2 = p.Results.Label2;
    figHandle = p.Results.FigureHandle;
    axHandle = p.Results.AxesHandle;
    
    % --- Create figure and axes ---
    if isempty(figHandle)
        figHandle = figure('Position', [100, 100, 800, 600], 'Name', 'Wing Planform Comparison');
    else
        figure(figHandle);
    end
    
    if isempty(axHandle)
        ax = axes('Parent', figHandle);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');
    else
        ax = axHandle;
        hold(ax, 'on');
    end
    
    % --- Get wing outlines for both designs ---
    [wing_x1, wing_y1, LE_x1, LE_y1, TE_x1, TE_y1] = getWingOutline(wingDesign1);
    [wing_x2, wing_y2, LE_x2, LE_y2, TE_x2, TE_y2] = getWingOutline(wingDesign2);
    
    % --- Plot wing 1 ---
    fill(ax, wing_x1, wing_y1, wingColor1, 'FaceAlpha', alphaVal, ...
         'EdgeColor', 'none', 'DisplayName', label1);
    
    % Plot leading edge for wing 1
    plot(ax, LE_x1, LE_y1, 'yo-', 'LineWidth', lineWidth, 'MarkerFaceColor', 'y', ...
         'DisplayName', [label1 ' Leading Edge']);
    
    % Plot trailing edge for wing 1
    plot(ax, TE_x1, TE_y1, 'ro-', 'LineWidth', lineWidth, 'MarkerFaceColor', 'r', ...
         'DisplayName', [label1 ' Trailing Edge']);
    
    % --- Plot wing 2 ---
    fill(ax, wing_x2, wing_y2, wingColor2, 'FaceAlpha', alphaVal, ...
         'EdgeColor', 'none', 'DisplayName', label2);
    
    % Plot leading edge for wing 2 with square markers
    plot(ax, LE_x2, LE_y2, 'ms-', 'LineWidth', lineWidth, 'MarkerFaceColor', 'm', ...
         'DisplayName', [label2 ' Leading Edge']);
    
    % Plot trailing edge for wing 2 with square markers
    plot(ax, TE_x2, TE_y2, 'rs-', 'LineWidth', lineWidth, 'MarkerFaceColor', 'r', ...
         'DisplayName', [label2 ' Trailing Edge']);
    
    % --- Set axis properties ---
    % Calculate bounds from both wings
    all_x = [wing_x1, wing_x2, LE_x1, LE_x2, TE_x1, TE_x2];
    all_y = [wing_y1, wing_y2, LE_y1, LE_y2, TE_y1, TE_y2];
    
    x_min = min(all_x) - 1;
    x_max = max(all_x) + 1;
    y_min = min(all_y) - 1;
    y_max = max(all_y) + 1;
    
    axis(ax, [x_min, x_max, y_min, y_max]);
    axis(ax, 'equal');
    
    % --- Add labels and title ---
    xlabel(ax, 'x (m)', 'FontSize', 11);
    ylabel(ax, 'y (m)', 'FontSize', 11);
    title(ax, 'Top-down view of wing planform comparison', 'FontSize', 12);
    
    % --- Add grid ---
    grid(ax, 'on');
    
    % --- Add legend ---
    legend(ax, 'Location', 'best', 'FontSize', 9);
    
    hold(ax, 'off');
    
    % --- Display summary in command window ---
    fprintf('\n=== Wing Planform Comparison ===\n');
    fprintf('%s:\n', label1);
    fprintf('  Span: %.2f m, Area: %.2f m², AR: %.2f\n', ...
            2*wingDesign1.y_tip, wingDesign1.S, wingDesign1.AR);
    fprintf('  Root chord: %.2f m, Tip chord: %.2f m\n', ...
            wingDesign1.c_root, wingDesign1.c_tip);
    fprintf('%s:\n', label2);
    fprintf('  Span: %.2f m, Area: %.2f m², AR: %.2f\n', ...
            2*wingDesign2.y_tip, wingDesign2.S, wingDesign2.AR);
    fprintf('  Root chord: %.2f m, Tip chord: %.2f m\n', ...
            wingDesign2.c_root, wingDesign2.c_tip);
    fprintf('================================\n\n');
end

function [wing_x, wing_y, LE_x, LE_y, TE_x, TE_y] = getWingOutline(wingDesign)
    % Helper function to get wing outline coordinates for top-down view
    
    % Extract coordinates
    x_root = wingDesign.x_root;
    x_kink = wingDesign.x_kink;
    x_tip  = wingDesign.x_tip;
    
    y_root = wingDesign.y_root;
    y_kink = wingDesign.y_kink;
    y_tip  = wingDesign.y_tip;
    
    c_root = wingDesign.c_root;
    c_kink = wingDesign.c_kink;
    c_tip  = wingDesign.c_tip;
    
    % Calculate trailing edge coordinates
    x_te_root = x_root + c_root;
    x_te_kink = x_kink + c_kink;
    x_te_tip  = x_tip  + c_tip;
    
    % Leading edge coordinates
    LE_x = [x_root, x_kink, x_tip];
    LE_y = [y_root, y_kink, y_tip];
    
    % Trailing edge coordinates
    TE_x = [x_te_root, x_te_kink, x_te_tip];
    TE_y = [y_root, y_kink, y_tip];
    
    % Combine for plotting the wing outline (closed polygon)
    wing_x = [LE_x, fliplr(TE_x), LE_x(1)];
    wing_y = [LE_y, fliplr(TE_y), LE_y(1)];
end

% Your existing code
dvec = DesignVector();
% dvec.LE_sweep = 20/180*pi;
wingDesign = WingDesign(dvec);

% initializer = load("fmincon_2026-01-05_17-29-48/initializer2026-01-05_17-24-20.mat").initializer;
% fminconresults = load("fmincon_2026-01-05_17-29-48/2026-01-05_19-13-35fmincon_results.mat");

initializer = load("fmincon_2026-01-07_11-55-38\initializer2026-01-07_11-51-48.mat").initializer;

fminconresults = load("fmincon_2026-01-07_11-55-38\2026-01-07_12-46-42fmincon_results.mat");
x_opt_normalized = fminconresults.x_opt;
x_opt = x_opt_normalized .* initializer.optimizer.x0;

dvec_opt = DesignVector();
dvec_opt = dvec_opt.fromVector(x_opt);
wingDesign_opt = WingDesign(dvec_opt);

% Plot comparison
plotWingPlanformComparison(initializer.optimizer.wingDesign, wingDesign_opt, ...
    'Label1', 'Initial Wing', ...
    'Label2', 'Optimized Wing', ...
    'Color1', [0.6 0.8 1], ...    % Light blue for initial
    'Color2', [0.6 1 0.6], ...    % Light green for optimized
    'Alpha', 0.5)
    