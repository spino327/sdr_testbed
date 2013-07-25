function [] = plotOPTA11(figure1, performance)
%PLOTOPTA Summary of this function goes here
%   Detailed explanation goes here

%     ml = [0.00171649, 5.00172, 10.0017, 15.0017, 20.0017, 25.0017, 30.0017];%, 35.0017];%, 40.0017, 45.0017];
%     mmse = [3.03444, 6.2145, 10.4289, 15.1449, 20.0495, 25.018, 30.0075];%, 35.0039];%, 40.0026, 45.0021];
%     ml = [1.2811, 3.1607, 9.2289, 14.4392];
    ml = [5.00172, 10.0017, 15.0017, 20.0017, 25.0017];%, 35.0017];%, 40.0017, 45.0017];
    mmse = [6.2145, 10.4289, 15.1449, 20.0495, 25.018];%, 35.0039];%, 40.0026, 45.0021];

    CSNR=5:5:25;    % expressed in dB
    WB = [1 1];%[1 1];%; 2 1];%; 3 1]%; 4 1; 3 2; 1 2; 1 3];%bandwidth ratio
    OPTA = zeros(size(WB, 1), length(CSNR));

    for i=1:size(WB, 1)
        for j=1:length(CSNR)
            %OPTA(i, j)=10.*log10((1 + 10^(0.1*CSNR(j)))^(1/2));
            OPTA(i, j)=10.*log10((1 + 10^(0.1*CSNR(j)))^(WB(i, 2)/WB(i, 1)));
        end
    end

    % Create axes
%     axes1 = axes('Parent',figure1,'YTick',[0 10 20 30 40 50],...
%         'YGrid','on',...
%         'XTick',[0 10 20 30 40 50],...
%         'XGrid','on');
%     box(axes1,'on');
%     hold(axes1,'all');

    grid on
    hold on;
    plot(CSNR, OPTA, 'k');

    % set(gca,'Ytick',[0:10:50]);
    xlabel('CSNR (dB)');
    ylabel(['SDR (dB)']);
    % title('Curvas OPTA');
    %legend('OPTA 1:1', 'OPTA 2:1', 'OPTA 3:1', 'OPTA 4:1', ...
    %    'OPTA 3:2', 'OPTA 1:2', 'OPTA 1:3',  'Location','SouthEast');

    plot(CSNR, ml, 'o', 'Color', 'r', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
    plot(CSNR, mmse, '+', 'Color', 'b', 'MarkerSize', 5);
    plot(performance(:, 2), performance(:, 3), 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 5);

    hold off;

    legend('OPTA 1:1', 'ML Decoding 1:1 (Simulation)', 'MMSE Decoding 1:1 (Simulation)', 'ML Decoding 1:1 (Testbed)')
%     legend('OPTA 2:1', 'ML Decoding 2:1 (Simulation)', 'ML Decoding 2:1 (Testbed)')

end

