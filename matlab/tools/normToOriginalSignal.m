function [ fix ] = normToOriginalSignal( x, max, min )
%UNTITLED2 Summary of this function goes here
%   [ fix ] = normToOriginalSignal( x, orginal )
%       in: - x = signal to normalize
%           - max = maximum value original signal
%           - min = minimum value original signal
%       out: - fix = signal x normalized

    fix = x;

    for i = 1:length(x)

        % max
        if x(i) > max
            fix(i) = max;
       
        % min
        elseif x(i) < min
            fix(i) = min;
        end
    end
end

