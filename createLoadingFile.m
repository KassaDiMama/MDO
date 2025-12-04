function createLoadingFile(Res,fileName)
    arguments
        Res
        fileName string
    end
    
    y = Res.Wing.Yst;
    cl = Res.Wing.cl;
    cm = Res.Wing.cm_c4;
    chord = Res.Wing.chord;
    L = cl*chord;
    M = cm*chord;
    
    fid = fopen(fileName+'.load', 'w');

    if fid == -1
        error('Cannot open boeing.init for writing.');
    end
    for i = 1:length(y)
        fprintf(fid, string(y)+" "+string(L)+" "+string(M)+"\n");
    end
    fclose(fid);
end