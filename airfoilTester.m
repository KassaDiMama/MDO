function plotDatCoords(filename)
    % PLOT DAT FILE WITH X-Y COORDINATES
    %
    % Input:
    %   filename - .dat file containing two numeric columns [x  y]
    %
    % Example:
    %   plotDatCoords('airfoil.dat')

    % Load the data (assumes no header, only numeric data)
    data = load(filename);

    % Ensure data has 2 columns
    if size(data,2) ~= 2
        error('The .dat file must contain exactly two numeric columns: x and y.');
    end

    x = data(:,1);
    y = data(:,2);

    % Plot
    % figure;
    plot(x, y, 'LineWidth', 2);
    axis equal;
    hold on; grid on;
    xlabel('x');
    ylabel('y');
    title(['Plot of ', filename]);
end

function plotCSTairfoil(N1, N2, AU, AL, styleU, styleL, labelU, labelL)
% plotCSTairfoil  Plots an airfoil using CST parameterization
%
% Optional styling and legend labels

    CST_order_U = length(AU) - 1;
    CST_order_L = length(AL) - 1;

    if CST_order_U ~= CST_order_L
        error('AU and AL must have the same CST order.');
    end

    CST_order = CST_order_U;
    t = linspace(0, 1, 400);

    yu = CSTcurve(t, AU, N1, N2, CST_order);
    yl = CSTcurve(t, AL, N1, N2, CST_order);

    plot(t, yu, styleU, 'LineWidth', 2, 'DisplayName', labelU);
    plot(t, yl, styleL, 'LineWidth', 2, 'DisplayName', labelL);
end

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
figure; hold on; grid on;

N1 = 0.5;
N2 = 1;

% ---- Original design ----
dvec = DesignVector();
wingDesign = WingDesign(dvec);

plotCSTairfoil( ...
    N1, N2, ...
    wingDesign.AU, wingDesign.AL, ...
    'r-', 'r--', ...
    'Original Upper', 'Original Lower');

% ---- Optimized design ----
initializer = load("initializer.mat").initializer;
fminconresults = load("fmincon_results.mat");

x_opt = fminconresults.x_opt .* initializer.optimizer.x0;

dvec = DesignVector();
dvec = dvec.fromVector(x_opt);
wingDesign_new = WingDesign(dvec);

plotCSTairfoil( ...
    N1, N2, ...
    wingDesign_new.AU, wingDesign_new.AL, ...
    'b-', 'b--', ...
    'Optimized Upper', 'Optimized Lower');

axis equal;
xlabel('x');
ylabel('y');
title('CST Airfoil Comparison');
legend('Location','best');
hold off;
