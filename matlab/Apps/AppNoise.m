%Copyright (c) 2012, Universidad Industrial de Santander, Colombia
%University of Delaware
%All rights reserved.
%@author: Sergio Pino
%@author: Henry Arguello
%Website: http://www.eecis.udel.edu/
%emails  : sergiop@udel.edu - henarfu@udel.edu
%Date   : Feb, 2012

classdef AppNoise
    %APPNOISE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p
        verbose = false
    end
    
    methods
        
        % constructor
        function self = AppNoise(p)
            self.p = p;
        end
        
        function self = createVectorOfSamples(self)
            
            files = dir(self.p);
           
            n = fix(10000);
            samples = zeros(1, n/2);
            samples2 = zeros(1, n/2);
            mrx = zeros(1, 2);
            vrx = zeros(1, 2);
            
            k = 1;
            for i = length(files):-1:1%length(files)-20
        
               if files(i).isdir == 0 && isempty(strfind(files(i).name, 'raw')) &&...
                    isempty(strfind(files(i).name, 'performance')) &&...
                    isempty(strfind(files(i).name, '.DS_Store'))

                    file = files(i).name;
                    fprintf('i: file = %s\n', file);
                    
                    fileData = strcat(self.p, file);
                    rx = real(read_complex_binary(fileData));
                    
%                     figure(1); plot(rx(1:200000));
                    
                    % find problematic peaks 
                    
                    for m = 1:1
                    
                        fixData = abs(rx);
                        meanData = mean(fixData);
                        varData = var(fixData);
                        [peaks, loc] = findpeaks(fixData, 'MINPEAKHEIGHT', 3*sqrt(varData) + 3*meanData);

                        if ~isempty(peaks)
                            
                        
                            % killing right part
                            j = 1;
                            exit = false;
                            while max(loc)+j <= length(rx) && exit ~= true
                                if fixData(max(loc)+j) < 3*sqrt(varData) + 3*meanData
                                    exit = true;
                                end
                                j = j + 1;
                            end

                            rx = rx(max(loc)+j+10:end);
                            
                        end
                    
                    end
                    
%                     figure(2); plot(rx);
                    
                    mrx(k) = mean(rx);
                    vrx(k) = var(rx);
                    
                    % U [a, b]; a + (b-a).*rand(100,1);
%                     pos = fix(1 + (length(rx)-1).*rand(1, n));
                    pos = fix(length(rx)/2-n/2):fix(length(rx)/2)-1;
%                     pos =(pos);
                    samples(k, :) = rx(pos);
                    
                    pos = fix(length(rx)/2):fix(length(rx)/2+n/2)-1;
                    samples2(k, :) = rx(pos);
                    
                    fprintf('i: file = %s, m = %d, v = %d\n', file, mrx(k), vrx(k));
                    fprintf('i: ms = %d, vs = %d\n', mean(samples(k, :)), var(samples(k, :)));
                    fprintf('i: ms2 = %d, vs2 = %d\n', mean(samples2(k, :)), var(samples2(k, :)));
                    
                    if self.verbose
                    
                        figure(1); plotspec(rx, 1);
                        figure(2); histfit(rx);
%                     
                        pause;
                    
                    end
                    
                    k = k + 1;

                    clear 'rx', 'fixData';
                    
               end
               
            end
            
            H = vartest2(reshape(samples, 1, []), reshape(samples2, 1, []))
            [p,STATS] = vartestn(samples')
            [p,STATS] = vartestn(samples2')
        end
        
        function performance = bpskSimSampledNoise(self, alpha, khist) 
            
            % setup
            files = dir(self.p);
            N = 20; % upsampling factor
            gd = 5;
            ro = 20;
            
            % signals
            ps = importdata(sprintf('data/rcosine/srrc_up%i_gd%i_ro%i.mat', N, gd, ro)); % pulse shape
            ps = ps/sqrt(sum((ps).^2));

            data = importdata('data/a_k/a_k1E4.mat');
            
            % results
            performance = zeros(1, 6); % Es/N0 snr BER BERThe Ne Net
            
            % BPSK
            modUtilsObj = ModUtils(ps, N, []);
            
            [sym] = modUtilsObj.bpskMod(data, 0, 0);
            [st, es] = modUtilsObj.pulseShaping(sym);
            nsamples = length(st);
            
            % running the test
            for pos = 1:length(khist)
                
                % getting the noise
                while 1
                    
                    i = fix(1 + (length(files)-1)*rand(1));
                    
                    if files(i).isdir == 0 && isempty(strfind(files(i).name, 'raw')) &&...
                            isempty(strfind(files(i).name, 'performance')) &&...
                            isempty(strfind(files(i).name, '.DS_Store'))
                        
                        file = files(i).name;
                        fprintf('i: file = %s\n', file);
                        
                        fileData = strcat(self.p, file);
                        nt = real(read_complex_binary(fileData));
                        
                        % find problematic peaks
                        
                        for m = 1:1
                            
                            fixData = abs(nt);
                            meanData = mean(fixData);
                            varData = var(fixData);
                            [peaks, loc] = findpeaks(fixData, 'MINPEAKHEIGHT', 3*sqrt(varData) + 3*meanData);
                            
                            if ~isempty(peaks)
                                
                                
                                % killing right part
                                j = 1;
                                exit = false;
                                while max(loc)+j <= length(nt) && exit ~= true
                                    if fixData(max(loc)+j) < 3*sqrt(varData) + 3*meanData
                                        exit = true;
                                    end
                                    j = j + 1;
                                end
                                
                                nt = nt(max(loc)+j+10:end);
                                
                            end
                            
                        end
                        
                        if length(nt) >= nsamples
                        
                            nt = nt(1:nsamples);
                            break;
                            
                        end
                            
                    end
                    
                end
                
                % model rt = alpha*k*st + n(t)
                rt = (alpha*khist(pos)).*st + nt';
                
                % dem
                rk = modUtilsObj.matchedFiltering(real(rt)); % output of the matched filter and sampling
                
                % decision
                recv = modUtilsObj.decision(rk);
                
                % Performance
                % CSNR
                csnr = (pow(rk) / pow(modUtilsObj.matchedFiltering(real(nt)))) - 1;
%                 fprintf('i: csnr: %f\n', 10*log10(csnr));
                csnr1 = 10*log10( (var(rk) + mean(rk)^2) / (var(modUtilsObj.matchedFiltering(real(nt))) + mean(modUtilsObj.matchedFiltering(real(nt)))^2) - 1 )
%                 csnrTheo = 10*log10( (var((k*m*alpha).*data) + mean((k*m*alpha).*data)^2) / (var(noiseThe) + mean(noiseThe)^2) )

                csnrH = alpha^2*khist(pos)^2*pow(data)/pow(modUtilsObj.matchedFiltering(real(nt)));

                fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));

                csnr = csnrH;


                % BER
                % matched filter
                [number_of_errors, bit_error_rate] = biterr(data, recv);
                fprintf('number of errors = %i, bit error rate = %f\n', number_of_errors, bit_error_rate);
                theoryBer = (1/2)*erfc(sqrt(2*csnr)/sqrt(2));%(1/2)*erfc(sqrt(csnr)/sqrt(2)); % theoretical ber
                fprintf('theory BER = %f\n', theoryBer);
                
                % Es/N0 snr BER BERThe Ne Net
                performance(pos, 1) = 10*log10(csnr);% Es_No
                performance(pos, 2) = csnr;
                performance(pos, 3) = bit_error_rate;
                performance(pos, 4) = theoryBer;
                performance(pos, 5) = number_of_errors;
                performance(pos, 6) = theoryBer*length(data);
                
            end
            
            % Plotting
            performance = sortrows(performance, -1);
            
            simBer = performance(:, 3); % simulated ber
            theoryBer = performance(:, 4); % theoretical ber
            
            figure;
            semilogy(performance(:, 1), theoryBer','b.-');
            hold on;
%             semilogy(performance(:, 1), simBer','x-');
            hold off;
            %     axis([-10  10^-5 0.5])
            grid on
            legend('theory', 'simulation');
            xlabel('Eb/No, dB');
            ylabel('Bit Error Rate');
            title('BER curve for BPSK modulation using real sampled noise');
            
        end
        
        
%         function 
%             
%         end
%         
        
    end
    
end

