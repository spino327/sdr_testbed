function [ rx ] = killingProblematicPeaks( rx, percent)
% find problematic peaks [ rx ] = killingProblematicPeaks( rx )
% kill the peaks bigger than the thresould = 3*sqrt(varData) + 3*meanData
% input:
%   - rx = signal, could be complex or real
%   - percent = percent of data to use in the algorithm

    if 0 > percent || 1 < percent
       disp('i: using 10% as default');
       percent = 0.1;
    end
    
    rxt = rx(1:floor(percent*length(rx)));

    fixData = abs(rxt);
    meanData = mean(fixData);
    varData = var(fixData);
    [peaks, loc] = findpeaks(fixData, 'MINPEAKHEIGHT', 3*sqrt(varData) + 3*meanData);

    if ~isempty(peaks)

        % killing right part
        j = 1;
        exit = false;
        while max(loc)+j <= length(rxt) && exit ~= true
            if fixData(max(loc)+j) < 3*sqrt(varData) + 3*meanData
                exit = true;
            end
            j = j + 1;
        end

        rx = rx(max(loc)+j+10:end);

    end

end

