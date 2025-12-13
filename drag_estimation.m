
function drag_wing = drag_estimation(Res, visc, show)

    if nargin < 3 || isempty(show)
        show = false;  % Default plotting = false
    end

    if visc == 1
        drag_wing = Res.CDwing;
        if show == true
            % Ensure column vectors and interpolate section Cd to wing Y locations
            wingY = Res.Wing.Yst(:);
            wingCdi = Res.Wing.cdi(:);
            secY = Res.Section.Y(:);
            secCd = Res.Section.Cd(:);

            % Use 'pchip' or 'linear'; 'extrap' prevents NaN outside range (or remove/exclude as needed)
            secCd_on_wing = interp1(secY, secCd, wingY, 'pchip', 'extrap');

            % Add both drag components (now same length)
            CD_total = wingCdi + secCd_on_wing;

            % Plot res.wing.cdi and res.section.cd on the same plot
            figure;
            hold on;
            plot(Res.Wing.Yst, Res.Wing.cdi, 'DisplayName', 'Wing CDi');
            plot(Res.Section.Y, Res.Section.Cd, 'DisplayName', 'Section CD');
            plot(wingY, CD_total, 'DisplayName', 'Total CD');
            xlabel('Spanwise Location');
            ylabel('Drag Coefficient');
            title('Drag Coefficients Comparison');
            legend show;
            grid on;
            hold off;
        end
    else
        drag_wing = Res.CDiwing;
        if show == true
            figure;
            hold on;
            plot(Res.Wing.Yst, Res.Wing.cdi, 'DisplayName', 'Wing CDi');
            xlabel('Spanwise Location');
            ylabel('Drag Coefficient');
            title('Drag Coefficients Distribution');
            legend show;
            grid on;
            hold off;
        end
    end
end