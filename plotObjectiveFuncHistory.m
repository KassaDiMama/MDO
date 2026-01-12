function plotOptimizationHistory(folderPath)
    % PLOTOPTIMIZATIONHISTORY Plot optimization history from iter_*.mat files
    %
    % Input:
    %   folderPath - Path to folder containing iter_*.mat files
    
    % If no folder specified, use current folder
    if nargin < 1
        folderPath = '.';
    end
    
    % Get list of iter_*.mat files
    fileList = dir(fullfile(folderPath, 'iter_*.mat'));
    
    
    % Initialize arrays
    iterations = zeros(length(fileList), 1);
    fvals = zeros(length(fileList), 1);
    
    % Extract data from each file
    for i = 1:length(fileList)
        try
            % Load the file
            data = load(fullfile(folderPath, fileList(i).name));
            
            % Extract iteration number from filename
            [~, filename, ~] = fileparts(fileList(i).name);
            iterations(i) = sscanf(filename, 'iter_%d');
            
            % Extract objective function value
            if isfield(data, 'optimValues') && isfield(data.optimValues, 'fval')
                fvals(i) = data.optimValues.fval;
            else
                error('optimValues.fval not found in file %s', fileList(i).name);
            end
            
        catch ME
            fprintf('Error loading file %s: %s\n', fileList(i).name, ME.message);
            iterations(i) = NaN;
            fvals(i) = NaN;
        end
    end
    
    % Remove any NaN entries
    validIdx = ~isnan(iterations) & ~isnan(fvals);
    iterations = iterations(validIdx);
    fvals = fvals(validIdx);
    
    if isempty(iterations)
        error('No valid data found in the specified files');
    end
    
    % Sort by iteration number (in case files were loaded out of order)
    [iterations, sortIdx] = sort(iterations);
    fvals = fvals(sortIdx);
    
    % Plot the optimization history
    figure('Position', [100, 100, 800, 600]);
    
    % Main plot
    subplot(1,1,1);
    plot(iterations, fvals, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', 'b');
    grid on;
    xlabel('Iteration Number', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Objective Function Value', 'FontSize', 12, 'FontWeight', 'bold');
    title('Optimization History', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add data labels for last few points
    if length(iterations) <= 10
        for i = 1:length(iterations)
            text(iterations(i), fvals(i), sprintf('  %.4g', fvals(i)), ...
                 'FontSize', 9, 'VerticalAlignment', 'bottom');
        end
    else
        % Only label first, last, and minimum
        text(iterations(1), fvals(1), sprintf('  Start: %.4g', fvals(1)), ...
             'FontSize', 9, 'VerticalAlignment', 'bottom');
        text(iterations(end), fvals(end), sprintf('  End: %.4g', fvals(end)), ...
             'FontSize', 9, 'VerticalAlignment', 'bottom');
        [minVal, minIdx] = min(fvals);
        text(iterations(minIdx), minVal, sprintf('  Min: %.4g', minVal), ...
             'FontSize', 9, 'VerticalAlignment', 'bottom');
    end
    
    
    
    % Display summary statistics
    fprintf('\n=== Optimization History Summary ===\n');
    fprintf('Number of iterations: %d\n', length(iterations));
    fprintf('Initial value: %.6g\n', fvals(1));
    fprintf('Final value: %.6g\n', fvals(end));
    fprintf('Total improvement: %.6g\n', fvals(1) - fvals(end));
    fprintf('Relative improvement: %.2f%%\n', ...
            100 * (fvals(1) - fvals(end)) / abs(fvals(1)));
end



plotOptimizationHistory("fmincon_2026-01-07_11-55-38")