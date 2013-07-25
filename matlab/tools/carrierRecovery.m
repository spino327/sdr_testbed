function [ ll ] = carrierRecovery( rlc )
%CARRIERRECOVERY Using a Large carrier signal for get the peak of the
%carrier
%   Detailed explanation goes here

    N = length(rlc);

%     figure(1); plotspec(rlc, 1)

    fftrlc = fft(rlc);
    [m, imax] = max(abs(fftrlc));

%     ssf = (0:N-1)/(N);
%     figure(2); hold on; plot(abs(fftrlc)); stem(imax, m); hold off;
% 
%     freqL = ssf(imax)
%     phaseL = angle(fftrlc(imax))

    ll = -((N-1)-imax);
%     rlcfix = rlc.*exp(j*2*pi*ll*(0:N-1)/N)';

%     fftrlc = fft(rlcfix);
%     [m, imax] = max(abs(fftrlc))
%     figure(3); hold on; plot(abs(fftrlc)); stem(imax, m); hold off;
%     figure(4); plotspec(rlcfix, 1)

end

