function stop = outputFcn(x, optimValues, state)
    stop = false;

    persistent logDir

    switch state
        case 'init'
            % Create timestamped folder
            tstr = datestr(now,'yyyy-mm-dd_HH-MM-SS');
            logDir = fullfile(pwd, ['fmincon_' tstr]);
            mkdir(logDir);

            % Save initial state
            save(fullfile(logDir,'init.mat'), ...
                'x', 'optimValues');

        case 'iter'
            % Save every iteration
            fname = sprintf('iter_%05d.mat', optimValues.iteration);
            save(fullfile(logDir,fname), ...
                'x', 'optimValues');

        case 'done'
            % Save final state
            save(fullfile(logDir,'final.mat'), ...
                'x', 'optimValues');

            % Clear persistent variable
            clear logDir
    end
end
