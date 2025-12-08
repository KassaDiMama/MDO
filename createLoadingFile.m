function createLoadingFile(Res,fileName,rho,V)
    arguments
        Res
        fileName string
        rho
        V
    end
    
    y = Res.Wing.Yst;
    cl = Res.Wing.cl;
    cm = Res.Wing.cm_c4;
    chord = Res.Wing.chord;
    L = cl.*chord.*y*0.5*rho*V^2;
    M = cm.*chord.*chord.*y*0.5*rho*V^2;
    
    fid = fopen(fileName+'.load', 'w');

    if fid == -1
        error('Cannot open boeing.init for writing.');
    end
    for i = 1:length(y)
        disp(string(y(i)) + " " + string(L(i)) + " " + string(M(i)));

        % write the line to file
        fprintf(fid, "%.6f %.6f %.6f\n", y(i), L(i), M(i));
        
    end
    fclose(fid);
end