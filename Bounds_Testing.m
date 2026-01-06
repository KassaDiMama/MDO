initializer = load("FINAL_CORRECT_fmincon_2025-12-24_18-11-55\initializer2025-12-24_18-05-34.mat").initializer;
obj = initializer;  % Assign the loaded initializer to the object

function y = CSTcurve(t, A, N1, N2, n)

                % Class function
                C = t.^N1 .* (1 - t).^N2;

                % Shape function 
                S = zeros(size(t));
                for i = 0:n
                    S = S + nchoosek(n, i) .* t.^i .* (1 - t).^(n - i) .* A(i + 1);
                end

                y = C .* S;
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


N1 = 0.5;
            N2 = 1;
            CST_order = length(obj.AU) - 1;

            % --- parametric domain ---
            ts = linspace(0, 1, 10000);  % resolution for plotting
            ratios = linspace(0,3,10000);
            % --- compute upper and lower surfaces ---
            for ratio_index = 1:length(ratios)
                ratio = ratios(ratio_index);
                yu = CSTcurve(ts, obj.AU*(1-ratio), N1, N2, CST_order);
                yl = CSTcurve(ts, obj.AL * (1+ratio), N1, N2, CST_order);

                mask_u = yu(2:end-1);
                mask_l = yl(2:end-1);
                res = mask_u > mask_l;
                if sum(res) < size(res,2)
                    % display(ratio)
                    fprintf('Intersection occurred at ratio of %f\n', ratio);
                    AU_lower_bound = (1 - ratios(ratio_index-1));
                    AL_upper_bound = (1 + ratios(ratio_index-1));
                    
                    fprintf('Lower bound: AU = %f, AL = %f\n', AU_lower_bound, AL_upper_bound);
                    
                    fprintf('Lower Ratio(i-1): %f, Upper Ratio (i): %f\n', ratios(ratio_index-1), ratios(ratio_index));
                    fprintf('Difference: %f\n', ratios(ratio_index) - ratios(ratio_index-1));
                    break;
                end
            end
            ratios = linspace(0,20,500);
            for ratio_index = 1:length(ratios)
                ratio = ratios(ratio_index);
                yu = CSTcurve(ts, obj.AU*(1+ratio), N1, N2, CST_order);
                yl = CSTcurve(ts, obj.AL * (1-ratio), N1, N2, CST_order);

                mask_u = yu(2:end-1);
                mask_l = yl(2:end-1);
                res = mask_u > mask_l;
                if max(yu+yl)>0.3
                    % display(ratio)
                    fprintf('Intersection occurred at ratio of %s\n', num2str(ratio));
                    % disp(max(yu+yl));
                    
                    AU_upper_bound = (1 + ratios(ratio_index-1));
                    AL_lower_bound = (1 - ratios(ratio_index-1));
                    break;
                end
            end
            % disp(AU_upper_bound);
            % disp(AL_lower_bound);
            % fprintf("AU: %f",1/max(obj.AU));
            % fprintf("AL: %f",-1/max(obj.AL))
        

% AUm(1:4) = obj.AU(1:4)*AU_upper_bound;
% AUm(4:end) = obj.AU(4:end)*(AU_upper_bound+0);
% ALm(1:4) = obj.AL(1:4)*(AL_lower_bound);
% ALm(4:5) = obj.AL(4:5)*(AL_lower_bound);
% ALm(end) = obj.AL(end)*(AL_lower_bound);
AUm = obj.AU*AU_lower_bound;
ALm = obj.AL*AL_lower_bound;
disp(ALm)
disp(AUm)

disp(AU_upper_bound);
disp(obj.AL*AL_lower_bound);
disp(AL_lower_bound);


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

plotCSTairfoil(...
    N1,N2,...
    AUm,ALm, ...
    'b-', 'b--', ...
    'Optimized Upper', 'Optimized Lower');

axis equal;
xlabel('x');
ylabel('y');
title('CST Airfoil Comparison');
legend('Location','best');
hold off;