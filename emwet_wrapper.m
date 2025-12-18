function emwet_wrapper(wingDesign,const,fileName, W_to, W_zf)
    arguments
        wingDesign WingDesign
        const   Const
        fileName
        W_to
        W_zf
    end
    fid = fopen(fileName+'.init', 'w');

    if fid == -1
        error('Cannot open boeing.init for writing.');
    end
    
    fprintf(fid, string(round(W_to))+" "+string(round(W_zf))+"\n");
    fprintf(fid, string(const.n_max)+"\n");
    fprintf(fid, "%.2f %.2f %d %d\n", wingDesign.S, wingDesign.b_half, wingDesign.number_of_platforms, wingDesign.number_of_airfoils);
    fprintf(fid, "%.2f %s\n", wingDesign.y_root, "b737a");
    fprintf(fid, "%.2f %s\n", wingDesign.y_kink, "b737a");
    fprintf(fid, "%.2f %s\n", wingDesign.y_tip, "b737a");
    
    fprintf(fid, "%.2f %.2f %.2f %.2f %.2f %.2f\n", ...
    wingDesign.c_root, wingDesign.x_root, wingDesign.y_root, wingDesign.z_root, wingDesign.front_spar_pos,wingDesign.rear_spar_pos);
    fprintf(fid, "%.2f %.2f %.2f %.2f %.2f %.2f | chord, x, y, z, front spar and rear spar position\n", ...
    wingDesign.c_kink, wingDesign.x_kink, wingDesign.y_kink, wingDesign.z_kink, wingDesign.front_spar_pos,wingDesign.rear_spar_pos);
    fprintf(fid, "%.2f %.2f %.2f %.2f %.2f %.2f\n", ...
    wingDesign.c_tip, wingDesign.x_tip, wingDesign.y_tip, wingDesign.z_tip, wingDesign.front_spar_pos,wingDesign.rear_spar_pos);

    fprintf(fid, "%.2f %.2f\n", wingDesign.start_tank, wingDesign.end_tank);
    fprintf(fid, "%d\n", wingDesign.engine_each_wing);
    fprintf(fid, "%.2f %d\n",wingDesign.engine_location, wingDesign.engine_weight);

    fprintf(fid, "%.5e %.2f %.5e %.5e\n", wingDesign.E_mod, wingDesign.density, wingDesign.tensile_yield, wingDesign.compressive_yield);
    fprintf(fid, "%.5e %.2f %.5e %.5e\n", wingDesign.E_mod, wingDesign.density, wingDesign.tensile_yield, wingDesign.compressive_yield);
    fprintf(fid, "%.5e %.2f %.5e %.5e\n", wingDesign.E_mod, wingDesign.density, wingDesign.tensile_yield, wingDesign.compressive_yield);
    
    fprintf(fid, "%.2f %.2f\n", wingDesign.efficiency_factor, wingDesign.rib_pitch);
    fprintf(fid, "%d\n", 1);

    fclose(fid);
end