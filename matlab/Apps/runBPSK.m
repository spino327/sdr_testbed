function [targetCSNR, alpha, khist, performance] = runBPSK()
%RUNBPSK Application for process BPSK received raw data from the USRP
%   

    % output
    targetCSNR = 0;
    alpha = 0;
    khist = zeros(1);
    
    %% params
    clc;
    genK = true;
    show = false;
%     showHist = true;
%     isRecording = false;
    N = 20; % upsampling factor
    gd = 5;
    ro = 20;
    f = 1;
    T = (1/f)*N;
    preamble = 500;
    pos_pream = 0;
    error = 0.001;
    path = '~/Dropbox/Research/data/rx/12-28-12-3/';
    nbins = 20;
    files = dir(path);

    % signals
    ps = importdata(sprintf('data/rcosine/srrc_up%i_gd%i_ro%i.mat', N, gd, ro)); % pulse shape
    ps = ps/sqrt(sum((ps).^2));
    ks = importdata('data/ks/ks101.mat');
    
    data = importdata('data/a_k/a_k1E4.mat');
    
    % results
    performance = zeros(1, 6); % Es/N0 snr BER BERThe Ne Net
    %%
    % processing
    pos = 1;
    for i = length(files):-1:1%length(files)-20
        
        if files(i).isdir == 0 && isempty(strfind(files(i).name, 'raw')) &&...
                isempty(strfind(files(i).name, 'performance')) &&...
                isempty(strfind(files(i).name, '.DS_Store'))
            
            file = files(i).name;
            fprintf('i: file = %s\n', strcat(path,file));
            
            % Here it begins
    
            %% parameters related to model rt = k*m*alpha*xt + nt
            % k
            init = strfind(file, 'at');
            final = strfind(file, '_G');
            k = str2num(file(init+2:final-1));%0.013656;%0.043185;%0.5;0.013656;0.024285;%0.004318;%0.007679;
            m = 1;
            
            khist(pos) = k;
            
            % signal name
            signame = file(1:init-1);
            
            %% sine
            sine = read_complex_binary(strcat(path, 'alpha/', signame, 'at1.0_G0.dat'));

            fileData = strcat(path, file);
            dt = read_complex_binary(fileData);

            % large carrier
            rlc = read_complex_binary(strcat(path, 'carrier/', signame, 'at1.0_G0.dat'));
            
            % expected nframes to transmitt
            lframe = 10000;
            exp_nframes = floor(length(data)/lframe);
            rec_frames = 0;
            if (lframe*exp_nframes < length(data))
                exp_nframes = exp_nframes+1;
            end

            % instance variables
            bpskObj = ModUtils(ps, N, ks);
            header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
            smRec = SmartReceiver(header);

            %% ALPHA estimate alpha, asume alpha cyclostationary

            % Estimate alpha from a sine wave
            %             [sineFix, cfix] = smRec.frequencySync(sine(floor(length(sine)/4)-floor(length(sine)/8):floor(length(sine)/4)+floor(length(sine)/8)));
            %             alpha = max(abs(real(sineFix)));
            sineFix = abs(sine(ceil(1/3*end):ceil(2/3*end)));
            alpha = mean(findpeaks(sineFix));
            fprintf('i: alpha = %f, k = %f\n', alpha, k);

        %     if show
        %         figure; hold on; plot(sineFix); hold off;
        %     end

            %% carrier recovery
            ll = carrierRecovery(rlc);
        
            %% what we received
            dtF = dt.*exp(1i*2*pi*ll*(0:length(dt)-1)/length(dt))';
            dtF = dt(10000:end);
            
            nframes = floor(length(dtF)/((lframe+length(ks)+preamble + 10)*N));
            if (((lframe+length(ks)+preamble + 10)*N)*nframes < length(dtF))
                nframes = nframes+1;
            end

            isFirst = true;
            fixPos = 0;

            recv = [];
            nt = [];

            for fr=1:nframes

                % receiving the portion of the message
                len = (lframe+length(ks)+preamble + 10)*N;

                if exp_nframes == 1
                    len2 = ((length(data)-lframe*rec_frames)+length(ks)+preamble + 10)*N;
                    rt = dtF(1 + (fr-1)*len + fixPos:1 + (fr-1)*len + fixPos + len2);
                elseif exp_nframes > 1
                    rt = dtF(1 + (fr-1)*len + fixPos:fr*len + fixPos);
                else
                    break;
                end

                % matched filtering
                rtf = bpskObj.mFilter(rt);
                [rtf, cfix] = smRec.frequencySync(rtf);
                rt = rt*exp(1i*cfix);
                fprintf('i: cumulative frequency fix = %f\n', cfix);

                %% rx

                % timing sync
                % finding the complete pilot, timing sync
                if exp_nframes == 1
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence(real(rt), (length(data)-lframe*rec_frames)*N);
                else
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence(real(rt), lframe*N);
                end
%                 result = length(headstart);

                if isempty(headstart) || length(headstart) > 10
                    continue;
                else

                    wasFound = false;
                    for j = length(headstart):-1:1
                        hs = headstart(j);

                        if (pow(real(rt))/pow(imag(rt))-1)>=0.1

                            if isFirst

                                fprintf('i: is first\n');
                                isFirst = false;
                                fixPos = hs - preamble*N;
                                rt = dtF(1 + (fr-1)*len + fixPos:fr*len + fixPos);
        %                         [rt, cfix] = smRec.frequencySync(rt);

                                % matched filtering
                                rtf = bpskObj.mFilter(rt);
                                [rtf, cfix] = smRec.frequencySync(rtf);
                                rt = rt*exp(1i*cfix);
                                fprintf('i: cumulative frequency fix = %f\n', cfix);

                                hs = preamble*N + 1;

                            end

                            if hs - preamble*N == -1
                                hs = hs + 1;
                            end

                            if hs >= preamble*N
                                fprintf('i: headstart found at sample %i', hs);
                                headstart = hs;
                                loc = loc(j);
                                wasFound = true;

                                exp_nframes = exp_nframes - 1;
                                rec_frames = rec_frames + 1;
                                break;

                            end
                        end
                    end

                    if wasFound == false || length(rt) < headstart+(length(lframe)+length(ks))*N;
                        disp('i: data not used due to ')
                        if isFirst == false
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                        end
                        continue;
                    end

                    if show
                        figure('Name', 'testBERPilotPlots - Pilot sequence');
                        subplot(2, 1, 1); hold on; plot(real(rtf)); stem(headstart, rt(headstart), 'r'); title('s(t)'); hold off;
                        subplot(2, 1, 2); hold on; plot(ycorr); stem(loc, ycorr(loc), 'g'); title('xcorr')
                    end

                end

                disp('i: finding marker sequences - OK')

                % setting starting sampling point
                init = headstart + 1;

                if show
                    figure('Name', 'testBERPilotPlots-RT');
                    hold on; plot(real(rtf), 'r'); plot(imag(rtf), 'b'); stem(init, rtf(init), 'sg'); title('r(t)'); hold off;
                end

                %Remove unwanted portions(first few samples till the peak value)
                noise = rt(headstart - preamble*N + 1:headstart);
                rt = rt(init:end);

                %% fixing phase
                if ycorr(loc) < 0
                    disp('i: fixing phase by multiplaying by cos(pi)');
                    rt = rt*cos(pi);
                end

                % dem
                rk = bpskObj.matchedFiltering(real(rt)); % output of the matched filter and sampling
        %         rk = bpskObj.downSampling(real(rt));

                % header
                result = bpskObj.decision(rk(1:length(ks)));

        %         mlunit.assert_equals(ks, result', sprintf('the ks is not as expected, nError = %i', biterr(ks, result')));

                % message
                if exp_nframes == 0
                    result = rk(length(ks) + 1:length(ks) + 1 + (length(data)-lframe*(rec_frames-1))-1);
                else
                    result = rk(length(ks) + 1:length(ks) + 1 + lframe-1);
                end

                if show
                    figure('Name', 'frame of image');
                    imshow(reshape(result(1:floor(length(result)/256)*256), 256, []), []); title('st');
                end

                recv = [recv result'];
                nt = [nt noise'];

            end

            %% SDR Real test
            fprintf('\ni: SDR - Testing the output MF and alpha*x + n\n');

            % as the theory said
            vn = var(real(nt));
            mn = mean(real(nt));
            noiseThe = mn + sqrt(vn)*randn(size(data));

            fprintf('i: vn = %d, mn = %d\n', vn, mn);

        %     rkDef = (k*m*alpha).*data + noiseThe;
        %     rkDef = rkDef(1:length(rk));
        % 
        %     if showHist
        %         figure;
        %         subplot(2, 1, 1); histfit(rkDef, nbins); title('SDR Real test, using the model rt = alpha*x + n');
        %         subplot(2, 1, 2); histfit(rk, nbins); title('SDR Real test, real match filtering - rk');
        %     end

            % Performance
            % CSNR
            csnr = (pow(recv) / pow(bpskObj.matchedFiltering(real(nt)))) - 1;
            
            csnr1 = 10*log10( (var(rk) + mean(rk)^2) / (var(bpskObj.matchedFiltering(real(nt))) + mean(bpskObj.matchedFiltering(real(nt)))^2) - 1 )
            csnrTheo = 10*log10( (var((k*m*alpha).*data) + mean((k*m*alpha).*data)^2) / (var(noiseThe) + mean(noiseThe)^2) )
            csnrH = alpha^2*k^2*m^2*pow(data)/pow(bpskObj.matchedFiltering(real(nt)));
            
            fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));
            
            csnr = csnrH;
            
            % BER
            % matched filter
            rk = bpskObj.decision(recv);
            [number_of_errors, bit_error_rate] = biterr(data, rk);
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
            pos = pos + 1;
        % 
        %     figure1 = figure('Name', 'OPTA');
        %     clf(figure1);
        %     plotOPTA(figure1, performance);

        if genK
            targetCSNR(:,1) = 0:2:50;
            targetCSNR(1,1) = 1E-10;
            for k=1:length(targetCSNR)
    
                targetCSNR(k,2) = generateK(targetCSNR(k,1), csnr);%(10^csnr)^(1/10));
            end
        end
    % End

        end
    end
    
    % Plotting
    performance = sortrows(performance, -1);
    
    simBer = performance(:, 3); % simulated ber
    theoryBer = performance(:, 4); % theoretical ber
    
    figure;
    semilogy(performance(:, 1), theoryBer','b.-');
    hold on;
    semilogy(performance(:, 1), simBer','mx-');
    hold off;
%     axis([-10  10^-5 0.5])
    grid on
    legend('theory', 'real transmission');
    xlabel('Eb/No, dB');
    ylabel('Bit Error Rate');
    title('testBERPilotPlots - Bit error probability curve for BPSK modulation');
    
    clear 'data' 'dt' 'dtF' 'fxml' 'fxmlThe' 'header' 'xml' 'xmlThe' 'ycorr'
    clear 'ks' 'noise' 'noiseThe' 'nt' 'ps' 'recv' 'result' 'rk' 'rkDef' 'rt' 'sine' 'sineFix'

end

