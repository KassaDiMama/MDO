function [CDi_total, CDv_total] = Drag_coeff_from_spanwise(Res, S_ref)
% Interpolates section data to wing stations and integrates
    
    % Extract data
    y_wing = Res.Wing.Yst(:);
    cdi = Res.Wing.cdi(:);
    chord = Res.Wing.chord(:);
    y_section = Res.Section.Y(:);
    Cd_visc = Res.Section.Cd(:);
    
    % Interpolate viscous drag to wing stations
    Cd_visc_interp = interp1(y_section, Cd_visc, y_wing, 'pchip', 'extrap');
    
    % Integrate
    CDi_total = (1/S_ref) * trapz(y_wing, chord .* cdi);
    CDv_total = (1/S_ref) * trapz(y_wing, chord .* Cd_visc_interp);
end