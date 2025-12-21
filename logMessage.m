function logMessage(msg, logFile)
    %LOGMESSAGE Write messages to a log file efficiently
    % msg can be a string, char, or string array
    % logFile is the filename (default: 'output.log')
    disp(msg);
    % if nargin < 2 || isempty(logFile)
    %     logFile = 'output.log';
    % end

    % % Convert to string array if necessary
    % if ischar(msg) || iscell(msg)
    %     msg = string(msg);
    % end

    % % Ensure column string array
    % if isrow(msg)
    %     msg = msg(:);
    % end

    % % Open file in append mode
    % fid = fopen(logFile, 'a');
    % if fid == -1
    %     error('Cannot open log file: %s', logFile);
    % end

    % % Join all lines into a single char block and prepend timestamp
    % timestamp = string(datestr(now,'yyyy-mm-dd HH:MM:SS'));
    % fullText = strjoin(msg, newline);   % string
    % fullText = string(fullText);
    % if ismissing(fullText)
    %     fullText = "<missing log message>";
    % end

    % fprintf(fid, "[%s] %s\n", string(timestamp), fullText);

    % % disp(fullText);  % also display in command window
    % fprintf(fid, '[%s] %s\n', char(timestamp), char(fullText));  % convert to char

    % % Close file
    % fclose(fid);
end
