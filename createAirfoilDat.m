function [t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1, N2, AU,AL,fileName)
    function result = CST(t,A)
        cn = t^N1 * (1-t)^N2;
        s = 0;
        for i = 0:5
            si = nchoosek(5,0)*t^i*(1-t)^(5-i)*A(i+1);
            s=s+si;
        end
        result = cn*s;
    end
    points_per_side = 46;

    t_upper = [];
    y_upper = [];
    for index = 0:points_per_side
        t = index/points_per_side;
        y = CST(t, AU);

        t_upper = [t_upper, t];
        y_upper = [y_upper, y];
    end
    
    t_upper = flip(t_upper);
    y_upper = flip(y_upper);

    t_lower = [];
    y_lower = [];
    for index = 0:points_per_side
        t = index/points_per_side;
        y = CST(t,AL);

        t_lower = [t_lower, t];
        y_lower = [y_lower, y];
    end
    % Calculate upper surface points using a similar approach
    
    
    fid = fopen(fileName+'.dat', 'w');

    if fid == -1
        error('Cannot open boeing.init for writing.');
    end
    for index = 1:length(t_upper)
        t = t_upper(index);
        y = y_upper(index);
        fprintf(fid, string(t)+" "+string(y)+"\n");
    end
    for index = 1:length(t_lower)
        t = t_lower(index);
        y = y_lower(index);
        fprintf(fid, string(t)+" "+string(y)+"\n");
    end
    fclose(fid);
end



