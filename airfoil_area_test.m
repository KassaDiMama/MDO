function [total_area, partial_area] = calculateAirfoilAreaPartial(datFilename, start_percent, end_percent)
    % Calculate airfoil area with option for partial chord area
    % Inputs:
    %   datFilename - path to .dat file
    %   start_percent - start percentage (0-100), default 0
    %   end_percent - end percentage (0-100), default 100
    % Outputs:
    %   total_area - total area from 0-100% chord
    %   partial_area - area from start_percent to end_percent
    
    % Set defaults if not provided
    if nargin < 2
        start_percent = 0;
    end
    if nargin < 3
        end_percent = 100;
    end
    
    % Validate input percentages
    if start_percent < 0 || start_percent > 100
        error('start_percent must be between 0 and 100');
    end
    if end_percent < 0 || end_percent > 100
        error('end_percent must be between 0 and 100');
    end
    if start_percent >= end_percent
        error('start_percent must be less than end_percent');
    end
    
    % Convert percentages to decimal fractions
    start_frac = start_percent / 100;
    end_frac = end_percent / 100;
    
    fprintf('Calculating area from %.1f%% to %.1f%% chord\n', start_percent, end_percent);
    
    % Load data using the simple method
    try
        data = load(datFilename);
    catch
        data = dlmread(datFilename);
    end
    
    if isempty(data)
        imported = importdata(datFilename);
        if isstruct(imported)
            data = imported.data;
        else
            data = imported;
        end
    end
    
    if isempty(data) || size(data, 2) < 2
        error('Could not load valid data from file.');
    end
    
    % Extract coordinates
    x = data(:, 1);
    y = data(:, 2);
    
    % Find minimum x (leading edge)
    [~, minIdx] = min(x);
    
    % Split data
    upper_x = x(1:minIdx);
    upper_y = y(1:minIdx);
    lower_x = x(minIdx:end);
    lower_y = y(minIdx:end);
    
    % Reverse upper surface
    upper_x = flipud(upper_x);
    upper_y = flipud(upper_y);
    
    % Sort
    [upper_x, idx] = sort(upper_x);
    upper_y = upper_y(idx);
    [lower_x, idx] = sort(lower_x);
    lower_y = lower_y(idx);
    
    % Create fine common grid for interpolation
    x_common = linspace(0, 1, 1000);
    
    % Interpolate both surfaces
    upper_y_interp = interp1(upper_x, upper_y, x_common, 'pchip', 'extrap');
    lower_y_interp = interp1(lower_x, lower_y, x_common, 'pchip', 'extrap');
    
    % Calculate total area (0-100%)
    total_area = trapz(x_common, upper_y_interp - lower_y_interp);
    
    % Calculate partial area
    % Find indices for the partial chord region
    idx_partial = x_common >= start_frac & x_common <= end_frac;
    x_partial = x_common(idx_partial);
    
    if isempty(x_partial)
        warning('No points found in the specified chord range. Adjusting boundaries.');
        
        % Find closest points
        [~, start_idx] = min(abs(x_common - start_frac));
        [~, end_idx] = min(abs(x_common - end_frac));
        
        if start_idx >= end_idx
            start_idx = max(1, end_idx - 10);
        end
        
        x_partial = x_common(start_idx:end_idx);
        upper_partial = upper_y_interp(start_idx:end_idx);
        lower_partial = lower_y_interp(start_idx:end_idx);
    else
        upper_partial = upper_y_interp(idx_partial);
        lower_partial = lower_y_interp(idx_partial);
    end
    
    % Calculate partial area
    partial_area = trapz(x_partial, upper_partial - lower_partial);
    
    % Display results
    fprintf('\n========================================\n');
    fprintf('Airfoil Area Analysis\n');
    fprintf('========================================\n');
    fprintf('File: %s\n', datFilename);
    fprintf('Total area (0-100%% chord): %.8f\n', total_area);
    fprintf('Partial area (%.1f%%-%.1f%% chord): %.8f\n', start_percent, end_percent, partial_area);
    fprintf('Partial/Total ratio: %.2f%%\n', (partial_area/total_area)*100);
    fprintf('========================================\n');
    
    % Plot results
    plotAirfoilWithPartialArea(x_common, upper_y_interp, lower_y_interp, ...
                               x_partial, upper_partial, lower_partial, ...
                               total_area, partial_area, start_percent, end_percent);
end

function plotAirfoilWithPartialArea(x_common, upper_y, lower_y, ...
                                    x_partial, upper_partial, lower_partial, ...
                                    total_area, partial_area, start_percent, end_percent)
    % Create figure
    figure('Position', [100, 100, 1200, 500], 'Name', 'Airfoil Area Analysis');
    
    % Plot 1: Full airfoil with highlighted section
    subplot(1, 3, 1);
    
    % Plot full airfoil
    fill([x_common, fliplr(x_common)], [upper_y, fliplr(lower_y)], ...
         [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Total area');
    hold on;
    
    % Highlight partial section
    fill([x_partial, fliplr(x_partial)], ...
         [upper_partial, fliplr(lower_partial)], ...
         [0.3, 0.8, 0.3], 'FaceAlpha', 0.7, 'EdgeColor', 'k', ...
         'LineWidth', 1.5, 'DisplayName', sprintf('%.0f-%.0f%% area', start_percent, end_percent));
    
    % Plot airfoil surfaces
    plot(x_common, upper_y, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Upper surface');
    plot(x_common, lower_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Lower surface');
    
    % Mark chord boundaries
    xline(start_percent/100, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, ...
          'DisplayName', sprintf('%.0f%% chord', start_percent));
    xline(end_percent/100, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, ...
          'DisplayName', sprintf('%.0f%% chord', end_percent));
    
    xlabel('x/c');
    ylabel('y/c');
    title(sprintf('Airfoil Cross-section\nTotal Area: %.6f', total_area));
    legend('Location', 'best');
    axis equal;
    grid on;
    
    % Plot 2: Thickness distribution with highlighted section
    subplot(1, 3, 2);
    
    thickness = upper_y - lower_y;
    thickness_partial = upper_partial - lower_partial;
    
    % Plot full thickness
    plot(x_common, thickness, 'k-', 'LineWidth', 2, 'DisplayName', 'Thickness');
    hold on;
    
    % Highlight partial section
    area_fill = fill([x_partial, fliplr(x_partial)], ...
                     [thickness_partial, zeros(size(thickness_partial))], ...
                     [1, 0.5, 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'r', ...
                     'LineWidth', 1.5, 'DisplayName', 'Partial area');
    
    % Mark chord boundaries
    xline(start_percent/100, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    xline(end_percent/100, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    
    xlabel('x/c');
    ylabel('Thickness');
    title(sprintf('Thickness Distribution\nPartial Area: %.6f', partial_area));
    legend('Location', 'best');
    grid on;
    
    % Plot 3: Data points comparison - FIXED VERSION
    subplot(1, 3, 3);
    
    % Find actual data points in the partial region
    bar_data = [total_area, partial_area];
    bar_labels = {'Total Area (0-100%)', sprintf('Partial Area (%d-%d%%)', start_percent, end_percent)};
    
    % Create bar chart with explicit color control
    bars = barh(bar_data);
    
    % Set bar colors correctly - FIXED HERE
    bars.FaceColor = 'flat';  % Enable per-bar coloring
    bars.CData(1,:) = [0.7, 0.7, 0.7];  % Gray for total
    bars.CData(2,:) = [0.3, 0.8, 0.3];  % Green for partial
    
    % Add value labels
    for i = 1:length(bar_data)
        text(bar_data(i) + max(bar_data)*0.01, i, ...
             sprintf('%.6f', bar_data(i)), ...
             'VerticalAlignment', 'middle', 'FontWeight', 'bold');
    end
    
    set(gca, 'YTickLabel', bar_labels);
    xlabel('Area');
    title('Area Comparison');
    grid on;
    
    % Add ratio text
    ratio = (partial_area/total_area)*100;
    annotation('textbox', [0.02, 0.02, 0.3, 0.05], ...
               'String', sprintf('Partial/Total = %.1f%%', ratio), ...
               'FitBoxToText', 'on', 'BackgroundColor', 'w', ...
               'FontSize', 10, 'FontWeight', 'bold');
end
% Example usage:
[total_area, partial_area] = calculateAirfoilAreaPartial('withcomb135.dat', 20, 60);