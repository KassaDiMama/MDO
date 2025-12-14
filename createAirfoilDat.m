function [t_upper,y_upper,t_lower, y_lower] = createAirfoilDat(N1, N2, AU,AL,fileName,ts,CST_order)
    % function result = CST(t,A)
    %     cn = t.^N1 .* (1-t).^N2;
    %     s = 0;
    %     for i = 0:5
    %         si = nchoosek(5,0)*t^i*(1-t)^(5-i)*A(i+1);
    %         s=s+si;
    %     end
    %     result = cn*s;
    % end
    function result = CST(t,A)
        cn = t.^N1 .* (1-t).^N2;
        s = 0;
        CST_order = length(A)-1;
        for i = 0:CST_order
            s = s + nchoosek(CST_order,i) * t.^i .* (1-t).^(CST_order-i) .* A(i+1);
        end
        result = cn .* s;
    end
    saving = false;
    if nargin < 6 || isempty(ts)
        points_per_side = 205;
        ts = linspace(0,1,points_per_side+1);
        saving = true;
    end
    
    
    t_upper = flip(ts);     
    t_upper = t_upper(1:end-1);
    y_upper = CST(t_upper, AU);
    
    % Lower surface
    t_lower = ts;
    y_lower = CST(t_lower, AL);
    % Calculate upper surface points using a similar approach
    
    if saving
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
end



