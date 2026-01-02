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

function plotWing3D(wingDesign, varargin)
    % plotWing3D - Plots a 3D isometric view of the wing geometry with CST airfoils
    %
    % Inputs:
    %   wingDesign - struct containing wing geometric parameters:
    %       x_root, y_root, z_root: root coordinates [m]
    %       x_kink, y_kink, z_kink: kink coordinates [m]
    %       x_tip, y_tip, z_tip: tip coordinates [m]
    %       c_root, c_kink, c_tip: chords [m]
    %       twist: 1x3 array of twist angles at root, kink, tip [deg] (positive = nose up)
    %       S: wing area [m²]
    %       AR: aspect ratio
    %       TR: taper ratio
    %       MAC: mean aerodynamic chord [m]
    %       AU: upper surface CST coefficients
    %       AL: lower surface CST coefficients
    %       incidence: incidence angle [deg] (positive = nose up)
    %
    % Optional Name-Value pairs:
    %   'Color' - Wing surface color [R G B] (default: [0.6 0.8 1])
    %   'Label' - Display label for legend (default: 'Wing')
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
    addParameter(p, 'Color', [0.6 0.8 1]);
    addParameter(p, 'Label', 'Wing');
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
    
    wingColor = p.Results.Color;
    wingLabel = p.Results.Label;
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
    twist_angles = wingDesign.twist; % degrees (positive = nose up)
    
    % Get incidence angle (positive = nose up)
    incidence = wingDesign.incidence; % degrees (positive = nose up)
    dihedral = wingDesign.dihedral;

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
    
    % --- Get other parameters ---
    S = wingDesign.S;
    AR = wingDesign.AR;
    TR = wingDesign.TR;
    MAC = wingDesign.MAC;
    
    % --- Create figure and axes ---
    if isempty(figHandle)
        figHandle = figure('Position', [100, 100, 1200, 800], 'Name', 'Wing 3D Geometry');
    else
        figure(figHandle);
    end
    
    if isempty(axHandle)
        ax = axes('Parent', figHandle);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');
        % Set equal scaling
        axis(ax, 'equal');
        daspect(ax, [1 1 1]);  % Force equal data aspect ratio
    else
        ax = axHandle;
        hold(ax, 'on');
    end
    
    % Set 3D view
    view(ax, viewAngle);
    
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
    sections_twist = twist_angles; % [root, kink, tip] twist in degrees (positive = nose up)
    
    % Trailing edge x positions
    sections_x_te = [x_te_root, x_te_kink, x_te_tip];
    
    n_vertices = length(t_full);
    n_sections = 3;
    
    X = zeros(n_vertices, n_sections);
    Y = zeros(n_vertices, n_sections);
    Z = zeros(n_vertices, n_sections);
    
    % Convert incidence to radians (positive = nose up)
    inc_rad = deg2rad(incidence);
    
    % --- Create airfoil sections with incidence AND twist ---
    for i = 1:n_sections
        x_le = sections_x_le(i);
        x_te = sections_x_te(i);
        y_section = sections_y(i);
        z_le = sections_z(i);
        chord = sections_chords(i);
        twist_deg = sections_twist(i); % Local twist at this section (positive = nose up)
        
        % Calculate total rotation angle for this section
        % Total rotation = incidence + local twist
        % Both are positive = nose up
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
                % Calculate vector from leading edge
                dx = x_unrotated - x_le;
                dz = z_unrotated - z_le;
                
                % Apply rotation matrix (positive angle = counterclockwise = nose up)
                X(j, i) = x_le + dx * cos(total_rotation_rad) - dz * sin(total_rotation_rad);
                Z(j, i) = z_le + dx * sin(total_rotation_rad) + dz * cos(total_rotation_rad);
            end
        end
        
        Y(:, i) = y_section * ones(n_vertices, 1);
    end
    
    % --- Plot the wing surface ---
    hWing = surf(ax, X, Y, Z, ...
        'FaceColor', wingColor, ...
        'FaceAlpha', alphaVal, ...
        'EdgeColor', edgeColor, ...
        'LineWidth', lineWidth, ...
        'DisplayName', wingLabel);
    
    % --- Plot airfoil sections if requested ---
    if showAirfoils
        colors = {'b', 'm', 'r'};
        labels = {'Root Airfoil', 'Kink Airfoil', 'Tip Airfoil'};
        
        for i = 1:n_sections
            plot3(ax, X(:, i), Y(:, i), Z(:, i), ...
                  [colors{i}, '-'], 'LineWidth', 2, 'DisplayName', labels{i});
        end
    end
    
    % --- Calculate actual LE and TE positions from airfoil surface ---
    % Find indices for LE and TE in the airfoil array
    le_idx = 1; % First point is leading edge
    [~, te_idx] = max(t_full); % Point with max t value (closest to 1)
    
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
    
    % --- Plot leading and trailing edges (using actual airfoil points) ---
    plot3(ax, LE_x_actual, LE_y_actual, LE_z_actual, ...
          'k-', 'LineWidth', 2.5, 'DisplayName', 'Leading Edge','HandleVisibility', 'off');
    plot3(ax, TE_x_actual, TE_y_actual, TE_z_actual, ...
          'k-', 'LineWidth', 2.5, 'DisplayName', 'Trailing Edge','HandleVisibility', 'off');
    
    % --- Plot chord lines at key stations ---
    station_colors = {'k', 'm', 'c'};
    station_names = {'Root', 'Kink', 'Tip'};
    chords = [c_root, c_kink, c_tip];
    
    for i = 1:3
        % Plot chord line (straight line between LE and TE)
        plot3(ax, [LE_x_actual(i), TE_x_actual(i)], ...
              [sections_y(i), sections_y(i)], ...
              [LE_z_actual(i), TE_z_actual(i)], ...
              [station_colors{i}, '--'], 'LineWidth', 1.5, ...
              'DisplayName', [station_names{i} ' Chord']);
        
        % Add chord labels (midpoint of chord line)
        mid_x = (LE_x_actual(i) + TE_x_actual(i)) / 2;
        mid_z = (LE_z_actual(i) + TE_z_actual(i)) / 2;
        
        label_text = sprintf('%s: %.2f m', station_names{i}, chords(i));
        if i == 2  % Special note for kink
            label_text = sprintf('%s: %.2f m', station_names{i}, chords(i));
        end
        
        text(ax, mid_x+1, sections_y(i), mid_z + 1, ...
             label_text, ...
             'HorizontalAlignment', 'center', 'FontSize', 9, ...
             'BackgroundColor', 'w', 'EdgeColor', station_colors{i});
    end
    
    % --- Add twist information to plot ---
    twist_text = sprintf('Twist: Root=%.1f°, Kink=%.1f°, Tip=%.1f°', ...
                         twist_angles(1), twist_angles(2), twist_angles(3));
    
    text(ax, min(X(:)) + (max(X(:))-min(X(:)))/2, ...
         max(Y(:)) + (max(Y(:))-min(Y(:)))*0.1, ...
         max(Z(:)) + (max(Z(:))-min(Z(:)))*0.1, ...
         twist_text, ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
         'BackgroundColor', 'w', 'EdgeColor', 'b');
    
    % --- Plot MAC if requested ---
    if showMAC
        % Calculate MAC spanwise position
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
        
        % Plot MAC line
        plot3(ax, [MAC_x_le, MAC_x_te], [y_bar, y_bar], [MAC_z_le, MAC_z_te], ...
              'g-', 'LineWidth', 2.5, 'DisplayName', 'MAC');
        
        % Plot MAC center point
        MAC_center_x = (MAC_x_le + MAC_x_te) / 2;
        MAC_center_z = (MAC_z_le + MAC_z_te) / 2;
        
        scatter3(ax, MAC_center_x, y_bar, MAC_center_z, 100, 'g', 'filled', ...
                 'DisplayName', 'MAC Center','HandleVisibility', 'off', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        
        % Label MAC
        text(ax, MAC_center_x+2, y_bar+0.5, MAC_center_z + 1, ...
             sprintf('MAC = %.2f m\nTwist = %.1f°', MAC, MAC_twist), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
             'BackgroundColor', 'w', 'EdgeColor', 'g');
    end
    
    % --- Plot key points ---
    % Leading edge points
    scatter3(ax, LE_x_actual, LE_y_actual, LE_z_actual, ...
             100, 'b', '^', 'filled', ...
             'DisplayName', 'Leading Edge Points', ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1);
    
    % Trailing edge points
    scatter3(ax, TE_x_actual, TE_y_actual, TE_z_actual, ...
             100, 'r', 'v', 'filled', ...
             'DisplayName', 'Trailing Edge Points', ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1);
    
    % --- Add labels to key points ---
    text(ax, LE_x_actual(1)-0.3, y_root - 0.3, LE_z_actual(1)+0.5, 'ROOT', ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
    text(ax, LE_x_actual(2)-0.3, y_kink - 0.3, LE_z_actual(2)+0.5, 'KINK', ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
    text(ax, LE_x_actual(3)-0.3, y_tip + 0.3, LE_z_actual(3)+0.5, 'TIP', ...
         'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
    
    % --- Plot coordinate system if requested ---
    if showCoordinates
        axis_length = max([y_tip, max([c_root, c_kink, c_tip])]) * 0.4;
        
        % X-axis
        plot3(ax, [0, axis_length], [0, 0], [0, 0], ...
              'k-', 'LineWidth', 2, 'DisplayName', 'X-axis','HandleVisibility', 'off');
        
        % Y-axis
        plot3(ax, [0, 0], [0, axis_length], [0, 0], ...
              'k-', 'LineWidth', 2, 'DisplayName', 'Y-axis','HandleVisibility', 'off');
        
        % Z-axis
        plot3(ax, [0, 0], [0, 0], [0, axis_length/3], ...
              'k-', 'LineWidth', 2, 'DisplayName', 'Z-axis','HandleVisibility', 'off');
        
        % Origin marker
        scatter3(ax, 0, 0, 0, 80, 'k', 'o', 'filled', ...
                 'DisplayName', 'Origin','HandleVisibility', 'off');
    end
    
    % --- Set axis properties ---
    % Calculate appropriate margins
    x_range = max(X(:)) - min(X(:));
    y_range = max(Y(:)) - min(Y(:));
    z_range = max(Z(:)) - min(Z(:));
    
    x_margin = max(x_range * 0.1, 0.5);
    y_margin = max(y_range * 0.1, 0.5);
    z_margin = max(z_range * 0.2, 0.2);
    
    x_min = min(X(:)) - x_margin;
    x_max = max(X(:)) + x_margin;
    y_min = min(Y(:)) - y_margin;
    y_max = max(Y(:)) + y_margin;
    z_min = min(Z(:)) - z_margin;
    z_max = max(Z(:)) + z_margin;
    
    % Ensure we include z=0 if near the wing
    if z_min > -0.1
        z_min = -0.2;
    end
    
    axis(ax, [x_min, x_max, y_min, y_max, z_min, z_max]);

    % Force equal scaling in all directions
    dataAspectRatio = [1 1 1];  % Change these if you want different scaling
    set(ax, 'DataAspectRatio', dataAspectRatio);
    set(ax, 'PlotBoxAspectRatio', [1 1 1]);  % Also set plot box aspect ratio
    
    % --- Add title and labels ---
    title_str = sprintf('Wing 3D Geometry');
    subtitle_str = sprintf('S = %.1f m² | AR = %.1f | dihedral = %.1f° | Incidence = %.1f°', ...
                          S, AR, dihedral, incidence);
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
    fprintf('\n=== Wing Geometry Summary ===\n');
    fprintf('Total span: %.2f m\n', 2*y_tip);
    fprintf('Area: %.2f m²\n', S);
    fprintf('Aspect Ratio: %.2f\n', AR);
    fprintf('Taper Ratio: %.3f\n', TR);
    fprintf('Mean Aerodynamic Chord: %.3f m\n', MAC);
    fprintf('Root chord: %.2f m\n', c_root);
    fprintf('Kink chord: %.2f m (TE aligned with root)\n', c_kink);
    fprintf('Tip chord: %.2f m\n', c_tip);
    fprintf('Incidence angle: %.1f° (positive = nose up)\n', incidence);
    fprintf('Twist: Root=%.1f°, Kink=%.1f°, Tip=%.1f° (positive = nose up)\n', twist_angles(1), twist_angles(2), twist_angles(3));
    fprintf('Dihedral at tip: %.1f°\n', atan2(z_tip - z_root, y_tip - y_root) * 180/pi);
    fprintf('CST Airfoils: Enabled (Order %d)\n', n);
    fprintf('=============================\n\n');
end



clear all
close all
clc


initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;

fminconresults = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\2025-12-24_20-15-52fmincon_results.mat");

x_opt_normalized = fminconresults.x_opt;
x_opt = x_opt_normalized .* initializer.optimizer.x0;

dvec = DesignVector();
dvec = dvec.fromVector(x_opt);
% dvec.LE_sweep = 20/180*pi;
wingDesign_new = WingDesign(dvec);

plotWing3D(wingDesign_new)