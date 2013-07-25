
function BPSKPlot(psim, preal)
    % BPSKPlot(psim, preal)
    % in: 
    %   - psim : simulation performance matrix [Es/N0, snr, BER, BERThe,
    %   Ne, Net]
    %   - preal : real performance matrix [Es/N0, snr, BER, BERThe,
    %   Ne, Net]
    % out: N/A

    % Plotting Performance
    psim = sortrows(psim, -1); % sim
    preal = sortrows(preal, -1); % real

    simBer = psim(:, 3); % simulated ber
    realBer = preal(:, 3); % real trasmission ber
    theoryBer = psim(:, 4); % theoretical ber

    figure;
    semilogy(psim(:, 1), theoryBer','b.-');
    hold on;
    semilogy(psim(:, 1), simBer','rx-');
    semilogy(preal(:, 1), realBer','ko-');
    hold off;
    %     axis([-10  10^-5 0.5])
    grid on
    legend('theory', 'simulation', 'real');
    xlabel('Eb/No, dB');
    ylabel('Bit Error Rate');
    title('BER curve for BPSK modulation');
    
end