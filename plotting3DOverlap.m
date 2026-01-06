function y = CSTcurve(t, A, N1, N2, n)
% CSTcurve  evaluates CST curve with given coefficients
%
% t  - parameter from 0 to 1
% A  - Bernstein coefficients
% N1, N2 - CST exponents
% n = CST order

    % Class function
    C = t.^N1 .* (1 - t).^N2;

    % Shape function (Bernstein polynomial)
    S = zeros(size(t));
    for i = 0:n
        S = S + nchoosek(n, i) .* t.^i .* (1 - t).^(n - i) .* A(i + 1);
    end

    y = C .* S;
end

function plotWing3D(wingDesign1, wingDesign2, varargin)
    % plotWing3D - Plots a 3D isometric view of two wing geometries with CST airfoils
    %
    % Inputs:
    %   wingDesign1 - struct containing first wing geometric parameters
    %   wingDesign2 - struct containing second wing geometric parameters
    %
    % Optional Name-Value pairs:
    %   'Color1' - First wing surface color [R G B] (default: [0.6 0.8 1])
    %   'Color2' - Second wing surface color [R G B] (default: [1 0.6 0.6])
    %   'Label1' - Display label for first wing legend (default: 'Initial Wing')
    %   'Label2' - Display label for second wing legend (default: 'Optimized Wing')
    %   'Alpha' - Face transparency (default: 0.6)
    %   'EdgeColor' - Edge color (default: 'none')
    %   'LineWidth' - Line width for edges (default: 1)
    %   'n_points' - Number of points for airfoil discretization (default: 50)
    %   'ShowAirfoils' - Show root and tip airfoil sections (default: true)
    %   'ShowMAC' - Show mean aerodynamic chord (default: true)
    %   'ShowCoordinates' - Show coordinate system (default: true)
    %   'ViewAngle' - View angle [azimuth elevation] (default: [45 30])
    %   'FigureHandle' - Handle to existing figure (default: creates new)
    %   'AxesHandle' - Handle to existing axes (default: creates new)
    
    % --- Defaults ---
    p = inputParser;
    addParameter(p, 'Color1', [0.6 0.8 1]);
    addParameter(p, 'Color2', [1 0.6 0.6]);
    addParameter(p, 'Label1', 'Initial Wing');
    addParameter(p, 'Label2', 'Optimized Wing');
    addParameter(p, 'Alpha', 0.6);
    addParameter(p, 'EdgeColor', 'none');
    addParameter(p, 'LineWidth', 1);
    addParameter(p, 'n_points', 50);
    addParameter(p, 'ShowAirfoils', true);
    addParameter(p, 'ShowMAC', true);
    addParameter(p, 'ShowCoordinates', true);
    addParameter(p, 'ViewAngle', [45 30]);
    addParameter(p, 'FigureHandle', []);
    addParameter(p, 'AxesHandle', []);
    parse(p, varargin{:});
    
    wingColor1 = p.Results.Color1;
    wingColor2 = p.Results.Color2;
    wingLabel1 = p.Results.Label1;
    wingLabel2 = p.Results.Label2;
    alphaVal = p.Results.Alpha;
    edgeColor = p.Results.EdgeColor;
    lineWidth = p.Results.LineWidth;
    n_points = p.Results.n_points;
    showAirfoils = p.Results.ShowAirfoils;
    showMAC = p.Results.ShowMAC;
    showCoordinates = p.Results.ShowCoordinates;
    viewAngle = p.Results.ViewAngle;
    figHandle = p.Results.FigureHandle;
    axHandle = p.Results.AxesHandle;
    
    % --- Create figure and axes ---
    if isempty(figHandle)
        figHandle = figure('Position', [100, 100, 1200, 800], 'Name', 'Wing 3D Geometry Comparison');
    else
        figure(figHandle);
    end
    
    if isempty(axHandle)
        ax = axes('Parent', figHandle);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');
        axis(ax, 'equal');
        daspect(ax, [1 1 1]);
    else
        ax = axHandle;
        hold(ax, 'on');
    end
    
    % Set 3D view
    view(ax, viewAngle);
    
    % --- Plot both wings ---
    % Get engine locations
    engine_y1 = wingDesign1.engine_location;
    engine_y2 = wingDesign2.engine_location;
    
    % Plot first wing (Initial Wing)
    plotSingleWing(ax, wingDesign1, wingColor1, wingLabel1, alphaVal, edgeColor, ...
                   lineWidth, n_points, showAirfoils, showMAC, false);
    
    % Plot second wing (Optimized Wing)
    plotSingleWing(ax, wingDesign2, wingColor2, wingLabel2, alphaVal, edgeColor, ...
                   lineWidth, n_points, showAirfoils, showMAC, false);
    
    % --- Plot engine locations for both wings ---
    % Find LE points for engine locations
    [LE_x1, LE_y1, LE_z1] = getWingPoints(wingDesign1, n_points);
    [LE_x2, LE_y2, LE_z2] = getWingPoints(wingDesign2, n_points);
    
    % Interpolate to find engine location on LE for first wing
    engine_x1 = interp1(LE_y1, LE_x1, engine_y1, 'linear', 'extrap');
    engine_z1 = interp1(LE_y1, LE_z1, engine_y1, 'linear', 'extrap');
    
    % Interpolate to find engine location on LE for second wing
    engine_x2 = interp1(LE_y2, LE_x2, engine_y2, 'linear', 'extrap');
    engine_z2 = interp1(LE_y2, LE_z2, engine_y2, 'linear', 'extrap');
    
    % Plot engine locations with different colors and symbols
    scatter3(ax, engine_x1, engine_y1, engine_z1, 150, 'r', 'pentagram', 'filled', ...
             'DisplayName', sprintf('%s Engine (y=%.1f m)', wingLabel1, engine_y1), ...
             'MarkerEdgeColor', 'm', 'LineWidth', 2);
    
    scatter3(ax, engine_x2, engine_y2, engine_z2, 150, 'b', 'hexagram', 'filled', ...
             'DisplayName', sprintf('%s Engine (y=%.1f m)', wingLabel2, engine_y2), ...
             'MarkerEdgeColor', 'y', 'LineWidth', 2);
    
    % Add text labels for engine locations with increased separation
    % Initial Wing engine label (positioned above and to the right)
    text(ax, engine_x1 -1.2, engine_y1 - 1, engine_z1 + 1.2, ...
         sprintf('%s Engine location\n%.1f m', wingLabel1, engine_y1), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'w', 'EdgeColor', 'm');
    
    % Optimized Wing engine label (positioned below and to the left)
    text(ax, engine_x2 - 1.2, engine_y2 + 1, engine_z2 + 1.5, ...
         sprintf('%s Engine location\n%.1f m', wingLabel2, engine_y2), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
         'BackgroundColor', 'w', 'EdgeColor', 'y');
    
    % --- Plot coordinate system if requested ---
    if showCoordinates
        % Calculate overall bounds from both wings
        [X1, Y1, Z1] = generateWingSurface(wingDesign1, n_points);
        [X2, Y2, Z2] = generateWingSurface(wingDesign2, n_points);
        
        X_all = [X1(:); X2(:)];
        Y_all = [Y1(:); Y2(:)];
        Z_all = [Z1(:); Z2(:)];
        
        axis_length = max([max(Y_all), max([wingDesign1.c_root, wingDesign1.c_kink, wingDesign1.c_tip, ...
                                            wingDesign2.c_root, wingDesign2.c_kink, wingDesign2.c_tip])]) * 0.4;
        
        % X-axis
        plot3(ax, [0, axis_length], [0, 0], [0, 0], ...
              'k-', 'LineWidth', 2, 'DisplayName', 'X-axis', 'HandleVisibility', 'off');
        
        % Y-axis
        plot3(ax, [0, 0], [0, axis_length], [0, 0], ...
              'k-', 'LineWidth', 2, 'DisplayName', 'Y-axis', 'HandleVisibility', 'off');
        
        % Z-axis
        plot3(ax, [0, 0], [0, 0], [0, axis_length/3], ...
              'k-', 'LineWidth', 2, 'DisplayName', 'Z-axis', 'HandleVisibility', 'off');
        
        % Origin marker
        scatter3(ax, 0, 0, 0, 80, 'k', 'o', 'filled', ...
                 'DisplayName', 'Origin', 'HandleVisibility', 'off');
    end
    
    % --- Set axis properties ---
    % Calculate appropriate margins from both wings
    [X1, Y1, Z1] = generateWingSurface(wingDesign1, n_points);
    [X2, Y2, Z2] = generateWingSurface(wingDesign2, n_points);
    
    X_all = [X1(:); X2(:)];
    Y_all = [Y1(:); Y2(:)];
    Z_all = [Z1(:); Z2(:)];
    
    x_range = max(X_all) - min(X_all);
    y_range = max(Y_all) - min(Y_all);
    z_range = max(Z_all) - min(Z_all);
    
    x_margin = max(x_range * 0.1, 0.5);
    y_margin = max(y_range * 0.1, 0.5);
    z_margin = max(z_range * 0.2, 0.2);
    
    x_min = min(X_all) - x_margin;
    x_max = max(X_all) + x_margin;
    y_min = min(Y_all) - y_margin;
    y_max = max(Y_all) + y_margin;
    z_min = min(Z_all) - z_margin;
    z_max = max(Z_all) + z_margin;
    
    % Ensure we include z=0 if near the wing
    if z_min > -0.1
        z_min = -0.2;
    end
    
    axis(ax, [x_min, x_max, y_min, y_max, z_min, z_max]);
    
    % Force equal scaling in all directions
    set(ax, 'DataAspectRatio', [1 1 1]);
    set(ax, 'PlotBoxAspectRatio', [1 1 1]);
    
    % --- Add title and labels ---
    title_str = 'Wing 3D Geometry Comparison';
    subtitle_str = sprintf('%s vs %s', wingLabel1, wingLabel2);
    title(ax, {title_str; subtitle_str}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(ax, 'X [m]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(ax, 'Y [m]', 'FontSize', 11, 'FontWeight', 'bold');
    zlabel(ax, 'Z [m]', 'FontSize', 11, 'FontWeight', 'bold');
    
    % --- Add grid ---
    grid(ax, 'on');
    grid(ax, 'minor');
    
    % --- Add legend ---
    legend(ax, 'Location', 'best', 'FontSize', 9);
    
    hold(ax, 'off');
    
    % Enable 3D rotation
    rotate3d(ax, 'on');
    
    % --- Display summary in command window ---
    fprintf('\n=== Wing Geometry Comparison ===\n');
    fprintf('%s:\n', wingLabel1);
    fprintf('  Span: %.2f m, Area: %.2f m², AR: %.2f\n', 2*wingDesign1.y_tip, wingDesign1.S, wingDesign1.AR);
    fprintf('  Root chord: %.2f m, Tip chord: %.2f m\n', wingDesign1.c_root, wingDesign1.c_tip);
    fprintf('  Incidence: %.1f°, Engine y=%.1f m\n', wingDesign1.incidence, engine_y1);
    fprintf('%s:\n', wingLabel2);
    fprintf('  Span: %.2f m, Area: %.2f m², AR: %.2f\n', 2*wingDesign2.y_tip, wingDesign2.S, wingDesign2.AR);
    fprintf('  Root chord: %.2f m, Tip chord: %.2f m\n', wingDesign2.c_root, wingDesign2.c_tip);
    fprintf('  Incidence: %.1f°, Engine y=%.1f m\n', wingDesign2.incidence, engine_y2);
    fprintf('================================\n\n');
end

function plotSingleWing(ax, wingDesign, wingColor, wingLabel, alphaVal, edgeColor, ...
                        lineWidth, n_points, showAirfoils, showMAC, showLegend)
    % Helper function to plot a single wing
    
    % --- Extract geometry ---
    x_root = wingDesign.x_root;
    x_kink = wingDesign.x_kink;
    x_tip  = wingDesign.x_tip;
    
    y_root = wingDesign.y_root;
    y_kink = wingDesign.y_kink;
    y_tip  = wingDesign.y_tip;

    z_root = wingDesign.z_root;
    z_kink = wingDesign.z_kink;
    z_tip  = wingDesign.z_tip;

    c_root = wingDesign.c_root;
    c_kink = wingDesign.c_kink;
    c_tip  = wingDesign.c_tip;

    % Get twist angles (1x3 array: [root, kink, tip])
    twist_angles = wingDesign.twist;
    
    % Get incidence angle
    incidence = wingDesign.incidence;
    
    % --- Trailing edge positions ---
    x_te_root = x_root + c_root;
    x_te_kink = x_root + c_root;  % Aligned with root TE
    x_te_tip = x_tip + c_tip;
    
    % --- CST parameters ---
    AU = wingDesign.AU;
    AL = wingDesign.AL;
    N1 = 0.5;
    N2 = 1.0;
    n = length(AU) - 1;
    
    % --- Generate CST airfoil sections with incidence and twist ---
    t = linspace(0, 1, n_points)';
    
    % Generate baseline airfoil (no incidence, no twist)
    y_upper = CSTcurve(t, AU, N1, N2, n);
    y_lower = CSTcurve(t, AL, N1, N2, n);
    
    % Create closed airfoil contour (no duplicate points)
    y_airfoil = [y_upper; flipud(y_lower(1:end-1))];
    t_full = [t; flipud(t(1:end-1))];
    
    % Define wing sections (root, kink, tip)
    sections_y = [y_root, y_kink, y_tip];
    sections_x_le = [x_root, x_kink, x_tip];
    sections_z = [z_root, z_kink, z_tip];
    sections_chords = [c_root, c_kink, c_tip];
    sections_twist = twist_angles;
    
    % Trailing edge x positions
    sections_x_te = [x_te_root, x_te_kink, x_te_tip];
    
    n_vertices = length(t_full);
    n_sections = 3;
    
    X = zeros(n_vertices, n_sections);
    Y = zeros(n_vertices, n_sections);
    Z = zeros(n_vertices, n_sections);
    
    % Convert incidence to radians
    inc_rad = deg2rad(incidence);
    
    % --- Create airfoil sections with incidence AND twist ---
    for i = 1:n_sections
        x_le = sections_x_le(i);
        x_te = sections_x_te(i);
        y_section = sections_y(i);
        z_le = sections_z(i);
        chord = sections_chords(i);
        twist_deg = sections_twist(i);
        
        % Calculate total rotation angle for this section
        total_rotation_deg = incidence + twist_deg;
        total_rotation_rad = deg2rad(total_rotation_deg);
        
        % For each point on the airfoil
        for j = 1:n_vertices
            % Normalized chord position
            chord_pos = t_full(j);
            
            % Original unrotated position
            x_unrotated = x_le + chord_pos * (x_te - x_le);
            z_unrotated = z_le + y_airfoil(j) * chord;
            
            if total_rotation_deg == 0
                % No rotation
                X(j, i) = x_unrotated;
                Z(j, i) = z_unrotated;
            else
                % Rotate around leading edge (x_le, z_le)
                dx = x_unrotated - x_le;
                dz = z_unrotated - z_le;
                
                % Apply rotation matrix
                X(j, i) = x_le + dx * cos(total_rotation_rad) - dz * sin(total_rotation_rad);
                Z(j, i) = z_le + dx * sin(total_rotation_rad) + dz * cos(total_rotation_rad);
            end
        end
        
        Y(:, i) = y_section * ones(n_vertices, 1);
    end
    
    % --- Plot the wing surface ---
    surf(ax, X, Y, Z, ...
        'FaceColor', wingColor, ...
        'FaceAlpha', alphaVal, ...
        'EdgeColor', edgeColor, ...
        'LineWidth', lineWidth, ...
        'DisplayName', wingLabel);
    
    % --- Plot airfoil sections if requested ---
    if showAirfoils
        % Determine colors based on wing type
        if contains(wingLabel, 'Initial')
            colors = {'r', 'm', 'c'};  % Red tones for initial wing
        else
            colors = {'b', 'g', 'y'};  % Blue tones for optimized wing
        end
        
        labels = {sprintf('%s Root', wingLabel), sprintf('%s Kink', wingLabel), sprintf('%s Tip', wingLabel)};
        
        for i = 1:n_sections
            plot3(ax, X(:, i), Y(:, i), Z(:, i), ...
                  [colors{i}, '-'], 'LineWidth', 2, 'DisplayName', labels{i});
        end
    end
    
    % --- Plot leading and trailing edges ---
    le_idx = 1;
    [~, te_idx] = max(t_full);
    
    LE_x_actual = zeros(1, 3);
    LE_y_actual = sections_y;
    LE_z_actual = zeros(1, 3);
    TE_x_actual = zeros(1, 3);
    TE_y_actual = sections_y;
    TE_z_actual = zeros(1, 3);
    
    for i = 1:3
        LE_x_actual(i) = X(le_idx, i);
        LE_z_actual(i) = Z(le_idx, i);
        TE_x_actual(i) = X(te_idx, i);
        TE_z_actual(i) = Z(te_idx, i);
    end
        
    % Plot edges without legend entries
    plot3(ax, LE_x_actual, LE_y_actual, LE_z_actual, ...
          'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot3(ax, TE_x_actual, TE_y_actual, TE_z_actual, ...
          'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    
    % --- Plot chord lines at key stations ---
    % Determine colors based on wing type
    if contains(wingLabel, 'Initial')
        station_colors = {'r', 'm', 'c'};  % Red tones for initial wing
    else
        station_colors = {'b', 'g', 'y'};  % Blue tones for optimized wing
    end
    
    station_names = {'Root', 'Kink', 'Tip'};
    chords = [c_root, c_kink, c_tip];
    
    for i = 1:3
        % Plot chord line (straight line between LE and TE)
        plot3(ax, [LE_x_actual(i), TE_x_actual(i)], ...
              [sections_y(i), sections_y(i)], ...
              [LE_z_actual(i), TE_z_actual(i)], ...
              [station_colors{i}, '--'], 'LineWidth', 1.5, ...
              'DisplayName', sprintf('%s %s Chord', wingLabel, station_names{i}),'HandleVisibility', 'off');
        
        % Add chord labels with offset based on wing type
        mid_x = (LE_x_actual(i) + TE_x_actual(i)) / 2;
        mid_z = (LE_z_actual(i) + TE_z_actual(i)) / 2;
        
        % Adjust label position based on wing type to separate them
        if contains(wingLabel, 'Initial')
            x_offset = 0.5;  % Right offset for initial wing
            z_offset = 1.5;  % Upward offset for initial wing
            y_offset = 0;
        else
            x_offset = 1.5; % Left offset for optimized wing
            z_offset = 1; % Downward offset for optimized wing
            y_offset = 0.8;
        end
        
        label_text = sprintf('%s: %.2f m', station_names{i}, chords(i));
        if i == 2  % Special note for kink
            label_text = sprintf('%s: %.2f m', station_names{i}, chords(i));
        end
        
        text(ax, mid_x + x_offset, sections_y(i)+ y_offset, mid_z + z_offset, ...
             label_text, ...
             'HorizontalAlignment', 'center', 'FontSize', 9, ...
             'BackgroundColor', 'w', 'EdgeColor', station_colors{i});
    end
    
    % --- Plot MAC if requested ---
    if showMAC
        % Calculate MAC spanwise position
        TR = wingDesign.TR;
        y_tip = wingDesign.y_tip;
        MAC = wingDesign.MAC;
        
        y_bar = y_tip * (1 + 2*TR) / (3*(1 + TR));
        
        % Interpolate to find MAC position and twist
        if y_bar <= y_kink
            % Inboard section
            alpha = (y_bar - y_root) / (y_kink - y_root);
            MAC_x_le = (1-alpha)*LE_x_actual(1) + alpha*LE_x_actual(2);
            MAC_z_le = (1-alpha)*LE_z_actual(1) + alpha*LE_z_actual(2);
            MAC_x_te = (1-alpha)*TE_x_actual(1) + alpha*TE_x_actual(2);
            MAC_z_te = (1-alpha)*TE_z_actual(1) + alpha*TE_z_actual(2);
            
            % Interpolate twist
            MAC_twist = (1-alpha)*twist_angles(1) + alpha*twist_angles(2);
        else
            % Outboard section
            alpha = (y_bar - y_kink) / (y_tip - y_kink);
            MAC_x_le = (1-alpha)*LE_x_actual(2) + alpha*LE_x_actual(3);
            MAC_z_le = (1-alpha)*LE_z_actual(2) + alpha*LE_z_actual(3);
            MAC_x_te = (1-alpha)*TE_x_actual(2) + alpha*TE_x_actual(3);
            MAC_z_te = (1-alpha)*TE_z_actual(2) + alpha*TE_z_actual(3);
            
            % Interpolate twist
            MAC_twist = (1-alpha)*twist_angles(2) + alpha*twist_angles(3);
        end
        
        % Use different MAC colors based on wing type
        if contains(wingLabel, 'Initial')
            mac_color = 'r';  % Red for initial wing
            mac_offset_x = 2.5;  % Offset to the right
            mac_offset_z = 1.5;  % Offset upward
        else
            mac_color = 'b';  % Blue for optimized wing
            mac_offset_x = -2.5; % Offset to the left
            mac_offset_z = -1.5; % Offset downward
        end
        
        % Plot MAC line
        plot3(ax, [MAC_x_le, MAC_x_te], [y_bar, y_bar], [MAC_z_le, MAC_z_te], ...
              [mac_color, '-'], 'LineWidth', 2.5, ...
              'DisplayName', sprintf('%s MAC (%.2f m)', wingLabel, MAC));
        
        % Plot MAC center point
        
        
    end
end

function [X, Y, Z] = generateWingSurface(wingDesign, n_points)
    % Helper function to generate wing surface points
    % Similar to plotSingleWing but returns points instead of plotting
    
    if nargin < 2
        n_points = 50;
    end
    
    % --- Extract geometry ---
    x_root = wingDesign.x_root;
    x_kink = wingDesign.x_kink;
    x_tip  = wingDesign.x_tip;
    
    y_root = wingDesign.y_root;
    y_kink = wingDesign.y_kink;
    y_tip  = wingDesign.y_tip;

    z_root = wingDesign.z_root;
    z_kink = wingDesign.z_kink;
    z_tip  = wingDesign.z_tip;

    c_root = wingDesign.c_root;
    c_kink = wingDesign.c_kink;
    c_tip  = wingDesign.c_tip;

    twist_angles = wingDesign.twist;
    incidence = wingDesign.incidence;
    
    % --- Trailing edge positions ---
    x_te_root = x_root + c_root;
    x_te_kink = x_root + c_root;
    x_te_tip = x_tip + c_tip;
    
    % --- CST parameters ---
    AU = wingDesign.AU;
    AL = wingDesign.AL;
    N1 = 0.5;
    N2 = 1.0;
    n = length(AU) - 1;
    
    % --- Generate CST airfoil sections ---
    t = linspace(0, 1, n_points)';
    
    y_upper = CSTcurve(t, AU, N1, N2, n);
    y_lower = CSTcurve(t, AL, N1, N2, n);
    y_airfoil = [y_upper; flipud(y_lower(1:end-1))];
    t_full = [t; flipud(t(1:end-1))];
    
    % Define wing sections
    sections_y = [y_root, y_kink, y_tip];
    sections_x_le = [x_root, x_kink, x_tip];
    sections_z = [z_root, z_kink, z_tip];
    sections_chords = [c_root, c_kink, c_tip];
    sections_twist = twist_angles;
    sections_x_te = [x_te_root, x_te_kink, x_te_tip];
    
    n_vertices = length(t_full);
    n_sections = 3;
    
    X = zeros(n_vertices, n_sections);
    Y = zeros(n_vertices, n_sections);
    Z = zeros(n_vertices, n_sections);
    
    inc_rad = deg2rad(incidence);
    
    % --- Generate points ---
    for i = 1:n_sections
        x_le = sections_x_le(i);
        x_te = sections_x_te(i);
        y_section = sections_y(i);
        z_le = sections_z(i);
        chord = sections_chords(i);
        twist_deg = sections_twist(i);
        
        total_rotation_deg = incidence + twist_deg;
        total_rotation_rad = deg2rad(total_rotation_deg);
        
        for j = 1:n_vertices
            chord_pos = t_full(j);
            x_unrotated = x_le + chord_pos * (x_te - x_le);
            z_unrotated = z_le + y_airfoil(j) * chord;
            
            if total_rotation_deg == 0
                X(j, i) = x_unrotated;
                Z(j, i) = z_unrotated;
            else
                dx = x_unrotated - x_le;
                dz = z_unrotated - z_le;
                X(j, i) = x_le + dx * cos(total_rotation_rad) - dz * sin(total_rotation_rad);
                Z(j, i) = z_le + dx * sin(total_rotation_rad) + dz * cos(total_rotation_rad);
            end
        end
        
        Y(:, i) = y_section * ones(n_vertices, 1);
    end
end

function [LE_x, LE_y, LE_z] = getWingPoints(wingDesign, n_points)
    % Helper function to get LE points for interpolation
    if nargin < 2
        n_points = 50;
    end
    
    % Generate airfoil points
    t = linspace(0, 1, n_points)';
    AU = wingDesign.AU;
    AL = wingDesign.AL;
    N1 = 0.5;
    N2 = 1.0;
    n = length(AU) - 1;
    
    y_upper = CSTcurve(t, AU, N1, N2, n);
    y_lower = CSTcurve(t, AL, N1, N2, n);
    y_airfoil = [y_upper; flipud(y_lower(1:end-1))];
    
    % Get sections
    sections_y = [wingDesign.y_root, wingDesign.y_kink, wingDesign.y_tip];
    sections_x_le = [wingDesign.x_root, wingDesign.x_kink, wingDesign.x_tip];
    sections_z = [wingDesign.z_root, wingDesign.z_kink, wingDesign.z_tip];
    sections_chords = [wingDesign.c_root, wingDesign.c_kink, wingDesign.c_tip];
    sections_twist = wingDesign.twist;
    
    x_te_root = wingDesign.x_root + wingDesign.c_root;
    x_te_kink = wingDesign.x_root + wingDesign.c_root;
    x_te_tip = wingDesign.x_tip + wingDesign.c_tip;
    sections_x_te = [x_te_root, x_te_kink, x_te_tip];
    
    LE_x = zeros(1, 3);
    LE_y = sections_y;
    LE_z = zeros(1, 3);
    
    inc_rad = deg2rad(wingDesign.incidence);
    
    for i = 1:3
        x_le = sections_x_le(i);
        z_le = sections_z(i);
        chord = sections_chords(i);
        twist_deg = sections_twist(i);
        
        total_rotation_deg = wingDesign.incidence + twist_deg;
        total_rotation_rad = deg2rad(total_rotation_deg);
        
        % LE point (t=0)
        x_unrotated = x_le;
        z_unrotated = z_le + y_airfoil(1) * chord;
        
        if total_rotation_deg ~= 0
            dx = x_unrotated - x_le;
            dz = z_unrotated - z_le;
            LE_x(i) = x_le + dx * cos(total_rotation_rad) - dz * sin(total_rotation_rad);
            LE_z(i) = z_le + dx * sin(total_rotation_rad) + dz * cos(total_rotation_rad);
        else
            LE_x(i) = x_unrotated;
            LE_z(i) = z_unrotated;
        end
    end
end

clear all
close all
clc


initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;

fminconresults = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\2025-12-24_20-15-52fmincon_results.mat");

% initializer = load("fmincon_2026-01-05_17-29-48\initializer2026-01-05_17-24-20.mat").initializer;
% 
% fminconresults = load("fmincon_2026-01-05_17-29-48\2026-01-05_19-13-35fmincon_results.mat");
x_opt_normalized = fminconresults.x_opt;
x_opt = x_opt_normalized .* initializer.optimizer.x0;

dvec = DesignVector();
dvec = dvec.fromVector(x_opt);
% dvec.LE_sweep = 20/180*pi;
wingDesign_new = WingDesign(dvec);

plotWing3D(initializer.optimizer.wingDesign,wingDesign_new)