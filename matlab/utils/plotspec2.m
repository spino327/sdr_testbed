% plotspec(x,Ts) plots the spectrum of the signal x
% Ts = time (in seconds) between adjacent samples in x

function plotspec2(x, Ts)

    N=length(x);                               % length of the signal x
    t=Ts*(1:N);                                % define a time vector 
    ssf=(-N/2:N/2-1)/(Ts*N);                   % frequency vector
    fx=fft(x(1:N));                            % do DFT/FFT
    fxs=fftshift(fx);                          % shift it for plotting
    subplot(3,1,1), hold on; plot(t, real(x), 'r'); plot(t, imag(x), 'b'); hold off;                % plot the waveform
    xlabel('seconds'); ylabel('amplitude')     % label the axes
    subplot(3,1,2), plot(ssf,abs(fxs))         % plot magnitude spectrum
    xlabel('frequency'); ylabel('magnitude')   % label the axes
    subplot(3,1,3), plot(ssf, angle(fxs));     % plot phase spectrum
    xlabel('frequency'); ylabel('angle')
    %verify parseval equalize power using
    %sum(abs(fxs).^2)/N
    %sum(abs(x).^2)
    % use axis([0,2,-1.1,1.1]) for specsquare.eps

    
end
