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
            ratios = linspace(1,0,1000);
            % --- compute upper and lower surfaces ---
            for ratio_index = 1:length(ratios)
                ratio = ratios(ratio_index);
                yu = CSTcurve(ts, obj.AU*ratio, N1, N2, CST_order);
                yl = CSTcurve(ts, obj.AL * ratio, N1, N2, CST_order);

                y_u = yu(2:end-1);
                y_l = yl(2:end-1);
                % res = mask_u > mask_l;
                if max(yu-yl)<Const.t_c_min
                    % display(ratio)
                    fprintf('Minimum thickness reached %f\n', max(yu-yl));
                    AU_lower_bound = ratios(round(ratio_index-1));
                    AL_lower_bound = ratios(round(ratio_index-1));
                    
                    fprintf('Lower bound: AU = %f, AL = %f\n', AU_lower_bound, AL_lower_bound);
                    mask_u = yu(2:end-1);
                    mask_l = yl(2:end-1);
                    res = mask_u < mask_l;
                    upperLowerOverlapFraction = sum(res)/length(mask_u);
                    
                    fprintf('Upper-Lower Overlap Fraction: %f\n', upperLowerOverlapFraction);
                    % fprintf('Lower Ratio(i-1): %f, Upper Ratio (i): %f\n', ratios(ratio_index-1), ratios(ratio_index));
                    fprintf('Difference: %f\n', ratios(ratio_index) - ratios(ratio_index-1));
                    break;
                end
            end
            ratios = linspace(1,5,1000);
            for ratio_index = 1:length(ratios)
                ratio = ratios(ratio_index);
                yu = CSTcurve(ts, obj.AU*ratio, N1, N2, CST_order);
                yl = CSTcurve(ts, obj.AL * ratio, N1, N2, CST_order);

                mask_u = yu(2:end-1);
                mask_l = yl(2:end-1);
                res = mask_u > mask_l;
                if max(yu-yl)>Const.t_c_max
                    fprintf('Maximum thickness reached %f\n', max(yu-yl));
                    AU_upper_bound = ratios(ratio_index-1);
                    AL_upper_bound = ratios(ratio_index-1);

                    mask_u = yu(2:end-1);
                    mask_l = yl(2:end-1);
                    res = mask_u < mask_l;
                    upperLowerOverlapFraction = sum(res)/length(mask_u);
                    
                    fprintf('Upper-Lower Overlap Fraction: %f\n', upperLowerOverlapFraction);
                    fprintf('Upper bound: AU = %f, AL = %f\n', AU_upper_bound, AL_upper_bound);
                    
                    % fprintf('Lower Ratio(i-1): %f, Upper Ratio (i): %f\n', ratios(ratio_index-1), ratios(ratio_index));
                    fprintf('Difference: %f\n', ratios(ratio_index) - ratios(ratio_index-1));
                    break
                end
            end
            % disp(AU_upper_bound);
            % disp(AU_lower_bound);
            % disp(AL_lower_bound);
            % disp(AL_upper_bound);
            % fprintf("AU: %f",1/max(obj.AU));
            % fprintf("AL: %f",-1/max(obj.AL))
        
      
AUm = obj.AU;
ALm = obj.AL;
disp(AUm);
i = 2;
AUm(i) = obj.AU(i)*AU_upper_bound;
disp(AUm);
% AUm(4:end) = obj.AU(4:end)*(AU_upper_bound+0);
% ALm(1:4) = obj.AL(1:4)*(AL_lower_bound);
% ALm(4:5) = obj.AL(4:5)*(AL_lower_bound);
% ALm(end) = obj.AL(end)*(AL_lower_bound);

% disp(ALm)
% disp(AUm)
% 
% disp(AU_upper_bound);
% disp(obj.AL*AL_lower_bound);
% disp(AL_lower_bound);


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