function [] = plotOPTA(figure1, performance, expML)
% PLOTOPTA Summary of this function goes here
% input:
%   - figure
%   - performance = matrix[n, 2] where 1st column is csnr and 2nd column is
%   sdr
%   - expML = matrix[n, 2] where 1st column is csnr and 2nd column is
%   sdr
    
%     performance = sortrows(performance, -1);
    CSNR = expML(:, 1);
%     expML = sortrows(expML, -1);
    
%     ml = [0.00171649, 5.00172, 10.0017, 15.0017, 20.0017, 25.0017, 30.0017];%, 35.0017];%, 40.0017, 45.0017];
%     mmse = [3.03444, 6.2145, 10.4289, 15.1449, 20.0495, 25.018, 30.0075];%, 35.0039];%, 40.0026, 45.0021];
%     ml = [1.2811, 3.1607, 9.2289, 14.4392];

%     CSNR=0:10:30;    % expressed in dB
    WB = [2 1];%[1 1];%; 2 1];%; 3 1]%; 4 1; 3 2; 1 2; 1 3];%bandwidth ratio
    OPTA = zeros(size(WB, 1), length(CSNR));

    for i=1:size(WB, 1)
        for j=1:length(CSNR)
            %OPTA(i, j)=10.*log10((1 + 10^(0.1*CSNR(j)))^(1/2));
            OPTA(i, j)=10.*log10((1 + 10^(0.1*CSNR(j)))^(WB(i, 2)/WB(i, 1)));
        end
    end

    % Create axes
    grid on
%     axes1 = axes('Parent',figure1,'YTick',[0:5:30],...
%         'YGrid','on',...
%         'XTick',[0:5:30],...
%         'XGrid','on');
%     box(axes1,'on');
%     hold(axes1,'all');

    hold on;
    plot(CSNR, OPTA, 'k');

    % set(gca,'Ytick',[0:10:50]);
    xlabel('CSNR (dB)');
    ylabel(['SDR (dB)']);
    % title('Curvas OPTA');
    %legend('OPTA 1:1', 'OPTA 2:1', 'OPTA 3:1', 'OPTA 4:1', ...
    %    'OPTA 3:2', 'OPTA 1:2', 'OPTA 1:3',  'Location','SouthEast');
    
%     for i=2:size(expML, 2)
        plot(expML(:, 1), expML(:, 2), 'o', 'Color', 'r', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
%         plot(expML(:, 1), expML(:, 3), '+', 'Color', 'r');
%         plot(expML(:, 1), expML(:, 4), '-s', 'Color', 'b');
%     end
    
%     plot(CSNR, mmse, '-o', 'Color', 'r');
    plot(performance(:, 1), performance(:, 2), 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 5);

    hold off;
    
%     legend('OPTA 1:1', 'ML Decoding 1:1 (Simulation)', 'MMSE Decoding 1:1 (Simulation)', 'ML Decoding 1:1 (Testbed)')
%     legend('OPTA 2:1', 'ML Decoding 2:1 (Simulation)', 'ML Decoding 2:1 (Testbed)')
%     legend('OPTA 2:1', 'ML Decoding (Sim) Image and testbed noise', 'ML Decoding (Sim) Image and gaussian noise', 'ML Decoding (Sim) pure gaussian', 'ML Decoding 2:1 (Testbed)')
%     legend('OPTA 2:1', 'ML Decoding (Sim) Image and testbed noise', 'ML Decoding (Sim) Image and gaussian noise', 'ML Decoding (Sim) pure gaussian')
%     legend('OPTA 2:1', 'ML Decoding (Sim) Image and testbed noise', 'ML Decoding (Sim) Image and gaussian noise')
    legend('OPTA 2:1', 'ML Decoding 2:1 (Simulation)', 'ML Decoding 2:1 (Testbed)')
end

