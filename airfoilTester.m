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
    figure;
    plot(x, y, 'LineWidth', 2);
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title(['Plot of ', filename]);
end

function plotCSTairfoil(N1, N2, AU, AL)
% plotCSTairfoil  Plots an airfoil using CST parameterization
%
% Inputs:
%   N1  - leading edge exponent
%   N2  - trailing edge exponent
%   AU  - upper surface CST coefficients
%   AL  - lower surface CST coefficients
%
% CST order is automatically detected:
%   CST_order = length(AU) - 1    (same for AL)
%
% Example:
%   plotCSTairfoil(0.5, 1.0, [0.1 0.2 0.05], [0 -0.05 -0.01])

    % --- automatic CST order detection ---
    CST_order_U = length(AU) - 1;
    CST_order_L = length(AL) - 1;

    if CST_order_U ~= CST_order_L
        error('AU and AL must have the same CST order.');
    end
    CST_order = CST_order_U;

    % --- parametric domain ---
    t = linspace(0, 1, 400);  % resolution for plotting

    % --- compute upper and lower surfaces ---
    yu = CSTcurve(t, AU, N1, N2, CST_order);
    yl = CSTcurve(t, AL, N1, N2, CST_order);

    % --- plot ---
    figure; hold on; grid on;
    plot(t, yu, 'r-', 'LineWidth', 2);
    plot(t, yl, 'b-', 'LineWidth', 2);
    axis equal;

    xlabel('x');
    ylabel('y');
    title('CST Airfoil');
    legend('Upper Surface','Lower Surface');

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
dvec = DesignVector();
wingDesign = WingDesign(dvec);
N1 = 0.5;
N2 = 1;
plotCSTairfoil(N1,N2,wingDesign.AU,wingDesign.AL);
plotDatCoords('airfoil.dat')