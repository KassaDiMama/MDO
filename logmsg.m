function logmsg(msg)
    disp(msg);
    return;
    persistent fid
    if isempty(fid)
        fname = "run_" + datestr(now,'yyyy-mm-dd_​HH-MM-SS') + ".txt";
        fid = fopen(fname,'a');
    end
    fprintf(fid, "%s\n", msg);
end