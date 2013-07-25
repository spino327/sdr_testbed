%Copyright (c) 2012, Universidad Industrial de Santander, Colombia
%University of Delaware
%All rights reserved.
%@author: Sergio Pino
%@author: Henry Arguello
%Website: http://www.eecis.udel.edu/
%emails  : sergiop@udel.edu - henarfu@udel.edu
%Date   : Feb, 2012

classdef AppBPSK
    %APPBPSK Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Instance variables
    properties
        N; % upsampling factor
        gd; % 
        ro; % roll of factor for the pulse shape filter
        sigpath; % path to signals
        data;
        n_preamble;
        pathrx; % path to the received signals
        ksf;
    end
    
    %% Instance methods
    methods
        
        % constructor
        function self = AppBPSK(N, gd, ro, sigpath, data, ksf, np, pathrx)
            self.N = N;
            self.gd = gd;
            self.ro = ro;
            self.sigpath = sigpath;
            self.data = data;
            self.n_preamble = np;
            self.pathrx = pathrx;
            self.ksf = ksf;
        end
        
        function st = createBSPKSignal(self, isLargeCarrier, signame, pathto)
            % creating a sample file of the modulated signal
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            a_k = importdata(strcat(self.sigpath, self.data));
            ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
            ks = (100/max(ps))*ones(1, length(ks));
            
            % mod
            bpskObj = ModUtils(ps, self.N, []);
            [sym] = [zeros(1, self.n_preamble) ks bpskObj.bpskMod(a_k, 0, 0)];
            st = bpskObj.pulseShaping(sym);
            
            % Praparing for small carrier or large carrier
            if(isLargeCarrier == 1)
                disp('Praparing large carrier modulation')
                st(self.n_preamble*self.N + 1:end) = st(self.n_preamble*self.N + 1:end) + 1;  
            end
            
%             save('data/carrier_recovery/bpskrlc.mat', 'st1', '-v7');
            if (isempty(pathto))
                write_float_binary(st, signame);
            else
                write_float_binary(st, strcat(pathto, signame));
            end
        end
            
%         function st = createBSPKSignal2(self, signame, pathto)
%             % creating a sample file of the modulated signal
%             
%             % signals
%             ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
%             ps = ps/sqrt(sum((ps).^2));
%             
%             a_k = importdata(strcat(self.sigpath, self.data));
%             ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
%             ks = (100/max(ps))*ones(1, length(ks));
%             
%             % mod
%             bpskObj = ModUtils(ps, self.N, []);
%             [sym] = [ks bpskObj.bpskMod(a_k, 0, 0)];
%             x = bpskObj.pulseShaping(sym);
%             
%             % m = |x + ashift|; ashift = |xmax - xmin|. ashift = amplitude
%             % shift
%             ashift = abs(max(x)-min(x));
%             
%             st = x + ashift;
%             
%             assert(sum(st) == sum(abs(st)));
%             
%             % Praparing for small carrier or large carrier
%             if(isLargeCarrier == 1)
%                 disp('Praparing large carrier modulation')
%                 st(self.n_preamble*self.N + 1:end) = st(self.n_preamble*self.N + 1:end) + 1;  
%             end
%             
% %             save('data/carrier_recovery/bpskrlc.mat', 'st1', '-v7');
%             if (isempty(pathto))
%                 write_float_binary(st, signame);
%             else
%                 write_float_binary(st, strcat(pathto, signame));
%             end
%         end
%         
        function [rk, initFix] = improveOutput(self, bpskObj, init, signal, data, factor, isLargeCarrier)
            
            if nargin < 4, factor = 0, isLargeCarrier = false; end
            
            pslen = floor(length(-self.gd:1/self.N:self.gd)/2);
            
            % left bound
            if init - pslen > 0, left = - pslen; else left = 1; end
            % right bound
            if init + pslen <= length(signal), right = pslen; else right = length(signal); end
            
            % [performance, sh]
            slideRes = zeros(length(left:right), 2);
            % Sliding from init - N to init + N, until get the best
            % decoding
            pos = 1;
            for sh = left:right
                
                sigTmp = signal(init + sh:end);
                
                % dem
                % output of the matched filter and sampling
                
                % removing large carrier
                if isLargeCarrier == true
                    sigTmp = bpskObj.matchedFiltering(abs(sigTmp) - factor);
                else
                    sigTmp = bpskObj.matchedFiltering(abs(sigTmp));
                end
                
                result = bpskObj.decision(sigTmp);
                
                if length(result) < length(data)
                    slideRes(pos, 1) = +Inf;
                    slideRes(pos, 2) = sh;
                    continue;
                end
                
                [number_of_errors, bit_error_rate, pos_errors] = biterr(data, result(1:length(data)));
                
                % [bit_error_rate, sh]
                slideRes(pos, 1) = bit_error_rate;
                slideRes(pos, 2) = sh;
                pos = pos + 1;
            end
            
            [v, index] = min(slideRes(:, 1));
            initFix = init + slideRes(index, 2);
            
            signal = signal(initFix : end);
            % removing large carrier
            if isLargeCarrier == true
                rk = bpskObj.matchedFiltering(abs(signal) - factor);
            else
                rk = bpskObj.matchedFiltering(abs(signal));
            end
        end
        
        function [target_csnr] = getKVector(self, csnr_vector, file, show)
            % getKVector return the vector of k values for the given CSNR
            % assuming a model y = alpha*k*x + n
            %
            
            % output
            
            % params
%             path = '~/Dropbox/Research/data/rx/12-28-12-3/';
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
            a_k = importdata(strcat(self.sigpath, self.data));

            % parameters related to model rt = k*fading*xt + nt
            
            % k
            init = strfind(file, 'at');
            final = strfind(file, '_G');
            k = str2double(file(init+2:final-1));
            
            % signal name
            signame = file(1:init-1);
            dt = read_complex_binary(strcat(self.pathrx, file));
            
            % sine
            sine = read_complex_binary(strcat(self.pathrx, 'fading/', signame, 'at1.0_G0.dat'));
            
            % expected nframes to receive
            lframe = 10000;
            exp_nframes = floor(length(a_k)/lframe);
            rec_frames = 0;
            
%             est_size = exp_nframes;
            
            if (lframe*exp_nframes < length(a_k))
                exp_nframes = exp_nframes + 1;
            end
            
            % instance variables
            bpskObj = ModUtils(ps, self.N, ks);
            header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
            smRec = SmartReceiver(header);
            
            % Fading estimate fading, asume fading cyclostationary
            fading = fadingEstimator(sine);
            
            fprintf('i: fading = %f, k = %f\n', fading, k);
            
            % carrier recovery
%             ll = carrierRecovery(dt);
            
            % what we received
%             dt = dt.*exp(1i*2*pi*ll*(0:length(dt)-1)/length(dt))';
            dt = dt(10000:end);
            
            nframes = floor( length(dt) / ((lframe+length(ks) + self.n_preamble + 10)*self.N) );
            if (((lframe + length(ks) + self.n_preamble + 10)*self.N ) * nframes < length(dt))
                nframes = nframes+1;
            end
            
            isFirst = true;
            fixPos = 0;
            
            recv = []; %zeros(1, est_size);
            nt = []; %zeros(1, est_size);
            
            for fr = 1:nframes
                
                % receiving the portion of the message
                len = (lframe + length(ks) + self.n_preamble + 10)*self.N;
                
                if exp_nframes == 1
                    len2 = ((length(a_k) - lframe*rec_frames) + length(ks) + self.n_preamble + 10)*self.N;
                    rt = dt(1 + (fr-1)*len + fixPos : 1 + (fr-1)*len + fixPos + len2);
                elseif exp_nframes > 1
                    rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                else
                    break;
                end
                
                % matched filtering
                rtf = bpskObj.mFilter(rt);
                [rtf, cfix] = smRec.frequencySync(rtf);
                rt = rt*exp(1i*cfix);
                fprintf('i: cumulative frequency fix = %f\n', cfix);
                
                % rx
                
                % timing sync
                % finding the complete pilot, timing sync
                if exp_nframes == 1
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence(real(rt), (length(a_k) - lframe*rec_frames)*self.N);
                else
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence(real(rt), lframe*self.N);
                end
                
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
                                fixPos = hs - self.n_preamble*self.N;
                                rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                                %                         [rt, cfix] = smRec.frequencySync(rt);
                                
                                % matched filtering
                                rtf = bpskObj.mFilter(rt);
                                [rtf, cfix] = smRec.frequencySync(rtf);
                                rt = rt*exp(1i*cfix);
                                fprintf('i: cumulative frequency fix = %f\n', cfix);
                                
                                hs = self.n_preamble*self.N + 1;
                                
                            end
                            
                            if hs - self.n_preamble*self.N == -1
                                hs = hs + 1;
                            end
                            
                            if hs >= self.n_preamble*self.N
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
                    
                    if wasFound == false || length(rt) < headstart+(length(lframe)+length(ks))*self.N;
                        disp('i: a_k not used due to ')
                        if isFirst == false
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                        end
                        continue;
                    end
                    
                    if show
                        figure('Name', 'r(t)');
                        subplot(2, 1, 1); hold on; plot(real(rt)); stem(headstart, rt(headstart), 'r'); title('r(t)'); hold off;
                        subplot(2, 1, 2); hold on; plot(ycorr); stem(loc, ycorr(loc), 'g'); title('xcorr')
                    end
                    
                end
                
                disp('i: finding marker sequences - OK')
                
                % setting starting sampling point
                init = headstart + 1;
                
                if show
                    figure('Name', 'testBERPilotPlots-RT');
                    hold on; plot(real(rt), 'r'); plot(imag(rt), 'b'); stem(init, rt(init), 'sg'); title('r(t)'); hold off;
                end
                
                %Remove unwanted portions(first few samples till the peak value)
                noise = rt(headstart - self.n_preamble*self.N + 1:headstart);
                rt = rt(init:end);
                
                %% fixing phase
                if ycorr(loc) < 0
                    disp('i: fixing phase by multiplaying by cos(pi)');
                    rt = rt*exp(1j*pi);
                end
                
                % dem
                % output of the matched filter and sampling
                rk = bpskObj.matchedFiltering(real(rt));
                
                % header
                rec_header = bpskObj.decision(rk(1:length(ks)));
                                
                % message
                if exp_nframes == 0
                    rec_sig = rk(length(ks) + 1 : length(ks) + 1 + (length(a_k)-lframe*(rec_frames-1))-1);
                else
                    rec_sig = rk(length(ks) + 1 : length(ks) + 1 + lframe-1);
                end
                
                recv = [recv rec_sig'];
                nt = [nt noise'];
                
            end
            
            %% SDR Real test
            fprintf('\ni: SDR - Testing the output MF and fading*x + n\n');
            
            % Ss the theory said
            vn = var(real(nt));
            mn = mean(real(nt));
            noiseThe = mn + sqrt(vn)*randn(size(a_k));
            
            fprintf('i: vn = %d, mn = %d\n', vn, mn);
            
            % Performance
            % CSNR
            csnr = (pow(recv) / pow(bpskObj.matchedFiltering(real(nt)))) - 1;
            
            csnr1 = 10*log10( (var(rk) + mean(rk)^2) / (var(bpskObj.matchedFiltering(real(nt))) + mean(bpskObj.matchedFiltering(real(nt)))^2) - 1 )
            csnrTheo = 10*log10( (var((k*fading).*a_k) + mean((k*fading).*a_k)^2) / (var(noiseThe) + mean(noiseThe)^2) )
            csnrH = fading^2*k^2*pow(a_k)/pow(bpskObj.matchedFiltering(real(nt)));
            
            fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));
            
            csnr = csnrH;
            
            target_csnr(:,1) = csnr_vector;
            for pos = 1 : length(target_csnr)
                
                target = target_csnr(pos, 1);
                if (target == 0)
                    target = 1E-16;
                end
                
                target_csnr(pos, 2) = k*generateK(target, csnr);
            end
            
        end
        
        function [target_csnr] = getKVector2(self, csnr_vector, file, isLargeCarrier, show)
            % getKVector return the vector of k values for the given CSNR
            % assuming a model y = alpha*k*x + n
            %
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
            ks = ones(1, length(ks));
            
            a_k = importdata(strcat(self.sigpath, self.data));

            % parameters related to model rt = k*fading*xt + nt
            
            % k
            init = strfind(file, 'at');
            final = strfind(file, '_G');
            k = str2double(file(init+2:final-1));
            
            % signal name
            signame = file(1:init-1);
            dt = read_complex_binary(strcat(self.pathrx, file));
            
            % sine
            sine = read_complex_binary(strcat(self.pathrx, 'fading/', signame, 'at1.0_G0.dat'));
            
            % expected nframes to receive
            lframe = 10000;
            exp_nframes = floor(length(a_k)/lframe);
            rec_frames = 0;
            
            if (lframe*exp_nframes < length(a_k))
                exp_nframes = exp_nframes + 1;
            end
            
            % instance variables
            bpskObj = ModUtils(ps, self.N, []);
            header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
            
            bpskObj = ModUtils(ps, self.N, ks);
            
            smRec = SmartReceiver(header);
            
            % Fading estimate fading, asume fading cyclostationary
            fading = fadingEstimator(sine);
            
            fprintf('i: fading = %f, k = %f\n', fading, k);
            
            if isLargeCarrier == true
                
                fixAmp = k*fading;
                header = header + 1;
            else
                
                fixAmp = 0;
               
            end
                
                
            dt = dt(10000:end);
            
            nframes = floor( length(dt) / ((lframe+length(ks) + self.n_preamble + 10)*self.N) );
            if (((lframe + length(ks) + self.n_preamble + 10)*self.N ) * nframes < length(dt))
                nframes = nframes+1;
            end
            
            isFirst = true;
            fixPos = 0;
            
            recv = [];
            nt = [];
            
            for fr = 1:nframes
                
                % receiving the portion of the message
                len = (lframe + length(ks) + self.n_preamble + self.gd*4)*self.N;
                
                if exp_nframes == 1
                    len2 = ((length(a_k) - lframe*rec_frames) + length(ks) + self.n_preamble + self.gd*4)*self.N;
                    rt = dt(1 + (fr-1)*len + fixPos : 1 + (fr-1)*len + fixPos + len2);
                elseif exp_nframes > 1
                    rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                else
                    break;
                end
                
                % rx
                
                % timing sync
                % finding the complete pilot, timing sync
                if exp_nframes == 1
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), (length(a_k) - lframe*rec_frames)*self.N);
                else
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), lframe*self.N);
                end
                
                if isempty(headstart) || length(headstart) > 10
                    continue;
                else
                    
                    wasFound = false;
                    for j = length(headstart):-1:1
                        hs = headstart(j);
                        
                        if isFirst == true
                            
                            fprintf('i: is first\n');
                            isFirst = false;
                            fixPos = hs - self.n_preamble*self.N;
%                             fixPos = hs - ((self.gd + 1) + self.n_preamble)*self.N;
                            rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos + 1);

                            hs = self.n_preamble*self.N + 1;
                            
                        end
                        
                        if hs - self.n_preamble*self.N == -1
                            hs = hs + 1;
                        end
                        
                        if hs >= self.n_preamble*self.N
                            fprintf('i: headstart found at sample %i', hs);
                            headstart = hs;
                            loc = loc(j);
                            wasFound = true;
                            
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                            break;
                            
                        end
                        
                    end
                    
                    if wasFound == false || length(rt) < headstart+(length(lframe)+length(ks))*self.N;
                        disp('i: a_k not used due to ')
                        if isFirst == false
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                        end
                        continue;
                    end
                    
                    if show == true
                        figure('Name', 'r(t)');
                        subplot(2, 1, 1); hold on; plot(real(rt)); stem(headstart, rt(headstart), 'r'); title('r(t)'); hold off;
                        subplot(2, 1, 2); hold on; plot(ycorr); stem(loc, ycorr(loc), 'g'); title('xcorr')
                    end
                    
                end
                
                disp('i: finding marker sequences - OK')

                %Remove unwanted portions(first few samples till the peak
                %value)
                % setting starting sampling point
%                 noise = rt(headstart - self.n_preamble*self.N + 1:headstart);
%                 rt = rt(init:end);
                
                initRt = headstart;
                initNoise = initRt - self.n_preamble*self.N;
                
                noise = rt(initNoise:initNoise + (self.n_preamble - 2*self.gd)*self.N);
                
                % larger carrier
%                 if isLargeCarrier == true
%                     initRt = initRt - (self.gd - 1)*self.N;
%                 end
                lenPS = length(-self.gd:1/self.N:self.gd);
                initRt = initRt - lenPS;
                                
                % fixing phase
                if ycorr(loc) < 0
                    disp('i: fixing phase by multiplaying by cos(pi)');
                    rt = rt*exp(1j*pi);
                end

                % Getting a more precise init of the Marker Sequence
                [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(...
                    abs(rt(initRt:initRt + 4*length(ks)*self.N)), 0);
                
                initRt = initRt + headstart - 1 - (self.gd-2)*self.N;
%                 rt = rt(initRt:end);
                
                % improving the performance of the frame
                if length(a_k) > lframe 
                    [rk, initRt] = self.improveOutput(bpskObj, initRt, rt, ...
                        [ks a_k((fr-1)*lframe + 1 : fr*lframe)]', fixAmp, isLargeCarrier);
                else
                    [rk, initRt] = self.improveOutput(bpskObj, initRt, rt, ...
                        [ks a_k]', fixAmp, isLargeCarrier);
                end

                
%                 %% fixing phase
%                 if ycorr(loc) < 0
%                     disp('i: fixing phase by multiplaying by cos(pi)');
%                     rt = rt*exp(1j*pi);
%                 end
%                 
%                 % dem
%                 % output of the matched filter and sampling
%                 
%                 % removing large carrier
%                 if isLargeCarrier == true
%                     rk = bpskObj.matchedFiltering(abs(rt) - k*fading);  
%                 else
%                     rk = bpskObj.matchedFiltering(abs(rt));
%                 end
                
                % header
                rec_header = bpskObj.decision(rk(1:length(ks)));
                                
                % message
                if exp_nframes == 0
                    rec_sig = rk(length(ks) + 1 : length(ks) + 1 + (length(a_k)-lframe*(rec_frames-1))-1);
                else
                    rec_sig = rk(length(ks) + 1 : length(ks) + 1 + lframe-1);
                end
                
                recv = [recv rec_sig'];
                nt = [nt noise'];
                
            end
            
            %% SDR Real test
            fprintf('\ni: SDR - Testing the output MF and fading*x + n\n');
            
            % Ss the theory said
            
            % Performance
            % CSNR
            csnr = (pow(recv) / pow(bpskObj.matchedFiltering(real(nt)))) - 1;
            csnrH = fading^2*k^2*pow(a_k)/pow(bpskObj.matchedFiltering(real(nt)));
%             csnr = (pow(recv) / pow(bpskObj.matchedFiltering(abs(nt) - k*fading))) - 1;
%             csnrH = fading^2*k^2*pow(a_k)/pow(bpskObj.matchedFiltering(abs(nt)-k*fading));
            fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));
            
            csnr = csnrH;
            
%             sk = bpskObj.bpskMod(a_k, 0, 0);
%             
            result = bpskObj.decision(recv);
%             
            [number_of_errors, bit_error_rate, pos_errors] = biterr(a_k, result);
            theoryBer = (1/2)*erfc(sqrt(2*csnr)/sqrt(2));
            fprintf('i: theory = %d, sim = %d, numerror = %i\n', theoryBer, bit_error_rate, number_of_errors);
%             %f, snr = %f, Es/No = %f', theoryBer, bit_error_rate, er_t, en, snr, Es_No));
            
            target_csnr(:,1) = csnr_vector;
            for pos = 1 : length(target_csnr)
                
                target = target_csnr(pos, 1);
                if (target == 0)
                    target = 1E-16;
                end
                
                target_csnr(pos, 2) = k*generateK(target, csnr);
            end
            
        end
        
        function [target_csnr] = getKVector3(self, csnr_vector, file, show)
            % getKVector return the vector of k values for the given CSNR
            % assuming a model y = alpha*k*x + n
            %
            
            % output
            
            % params
%             path = '~/Dropbox/Research/data/rx/12-28-12-3/';
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
            a_k = importdata(strcat(self.sigpath, self.data));

            % parameters related to model rt = k*fading*xt + nt
            
            % k
            init = strfind(file, 'at');
            final = strfind(file, '_G');
            k = str2double(file(init+2:final-1));
            
            % signal name
            signame = file(1:init-1);
            dt = read_complex_binary(strcat(self.pathrx, file));
            
            % sine
            sine = read_complex_binary(strcat(self.pathrx, 'fading/', signame, 'at1.0_G0.dat'));
            
            % expected nframes to receive
            lframe = 10000;
            exp_nframes = floor(length(a_k)/lframe);
            rec_frames = 0;
            
%             est_size = exp_nframes;
            
            if (lframe*exp_nframes < length(a_k))
                exp_nframes = exp_nframes + 1;
            end
            
            % instance variables
            bpskObj = ModUtils(ps, self.N, ks);
            header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
            smRec = SmartReceiver(header);
            
            % Fading estimate fading, asume fading cyclostationary
            fading = fadingEstimator(sine);
            
            fprintf('i: fading = %f, k = %f\n', fading, k);
            
            % carrier recovery
%             ll = carrierRecovery(dt);
            
            % what we received
%             dt = dt.*exp(1i*2*pi*ll*(0:length(dt)-1)/length(dt))';
            dt = dt(10000:end);
            
            nframes = floor( length(dt) / ((lframe+length(ks) + self.n_preamble + 10)*self.N) );
            if (((lframe + length(ks) + self.n_preamble + 10)*self.N ) * nframes < length(dt))
                nframes = nframes+1;
            end
            
            isFirst = true;
            fixPos = 0;
            
            recv = []; %zeros(1, est_size);
            nt = []; %zeros(1, est_size);
            
            for fr = 1:nframes
                
                % receiving the portion of the message
                len = (lframe + length(ks) + self.n_preamble + 10)*self.N;
                
                if exp_nframes == 1
                    len2 = ((length(a_k) - lframe*rec_frames) + length(ks) + self.n_preamble + 10)*self.N;
                    rt = dt(1 + (fr-1)*len + fixPos : 1 + (fr-1)*len + fixPos + len2);
                elseif exp_nframes > 1
                    rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                else
                    break;
                end
                
                % matched filtering
                rtf = bpskObj.mFilter(rt);
                [rtf, cfix] = smRec.frequencySync(rtf);
                rt = rt*exp(1i*cfix);
                fprintf('i: cumulative frequency fix = %f\n', cfix);
                
                % rx
                
                % timing sync
                % finding the complete pilot, timing sync
                if exp_nframes == 1
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), (length(a_k) - lframe*rec_frames)*self.N);
                else
                    [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), lframe*self.N);
                end
                
                if isempty(headstart) || length(headstart) > 10
                    continue;
                else
                    
                    wasFound = false;
                    for j = length(headstart):-1:1
                        hs = headstart(j);
                        
                        if isFirst
                            
                            fprintf('i: is first\n');
                            isFirst = false;
                            fixPos = hs - self.n_preamble*self.N;
                            rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                            %                         [rt, cfix] = smRec.frequencySync(rt);
                            
                            % matched filtering
                            rtf = bpskObj.mFilter(rt);
                            [rtf, cfix] = smRec.frequencySync(rtf);
                            rt = rt*exp(1i*cfix);
                            fprintf('i: cumulative frequency fix = %f\n', cfix);
                            
                            hs = self.n_preamble*self.N + 1;
                            
                        end
                        
                        if hs - self.n_preamble*self.N == -1
                            hs = hs + 1;
                        end
                        
                        if hs >= self.n_preamble*self.N
                            fprintf('i: headstart found at sample %i', hs);
                            headstart = hs;
                            loc = loc(j);
                            wasFound = true;
                            
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                            break;
                            
                        end
                        
                    end
                    
                    if wasFound == false || length(rt) < headstart+(length(lframe)+length(ks))*self.N;
                        disp('i: a_k not used due to ')
                        if isFirst == false
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                        end
                        continue;
                    end
                    
                    if show
                        figure('Name', 'r(t)');
                        subplot(2, 1, 1); hold on; plot(real(rt)); stem(headstart, rt(headstart), 'r'); title('r(t)'); hold off;
                        subplot(2, 1, 2); hold on; plot(ycorr); stem(loc, ycorr(loc), 'g'); title('xcorr')
                    end
                    
                end
                
                disp('i: finding marker sequences - OK')
                
                % setting starting sampling point
                init = headstart + 1;
                
                if show
                    figure('Name', 'testBERPilotPlots-RT');
                    hold on; plot(real(rt), 'r'); plot(imag(rt), 'b'); stem(init, rt(init), 'sg'); title('r(t)'); hold off;
                end
                
                %Remove unwanted portions(first few samples till the peak value)
                noise = rt(headstart - self.n_preamble*self.N + 1:headstart);
                rt = rt(init:end);
                
                %% fixing phase
                if ycorr(loc) < 0
                    disp('i: fixing phase by multiplaying by cos(pi)');
                    rt = rt*exp(1j*pi);
                end
                
                % dem
                % output of the matched filter and sampling
                rk = bpskObj.matchedFiltering(real(rt));
                
                % header
                rec_header = bpskObj.decision(rk(1:length(ks)));
                                
                % message
                if exp_nframes == 0
                    rec_sig = rk(length(ks) + 1 : length(ks) + 1 + (length(a_k)-lframe*(rec_frames-1))-1);
                else
                    rec_sig = rk(length(ks) + 1 : length(ks) + 1 + lframe-1);
                end
                
                recv = [recv rec_sig'];
                nt = [nt noise'];
                
            end
            
            %% SDR Real test
            fprintf('\ni: SDR - Testing the output MF and fading*x + n\n');
            
            % Ss the theory said
            vn = var(real(nt));
            mn = mean(real(nt));
            noiseThe = mn + sqrt(vn)*randn(size(a_k));
            
            fprintf('i: vn = %d, mn = %d\n', vn, mn);
            
            % Performance
            % CSNR
            csnr = (pow(recv) / pow(bpskObj.matchedFiltering(real(nt)))) - 1;
            
            csnr1 = 10*log10( (var(rk) + mean(rk)^2) / (var(bpskObj.matchedFiltering(real(nt))) + mean(bpskObj.matchedFiltering(real(nt)))^2) - 1 )
            csnrTheo = 10*log10( (var((k*fading).*a_k) + mean((k*fading).*a_k)^2) / (var(noiseThe) + mean(noiseThe)^2) )
            csnrH = fading^2*k^2*pow(a_k)/pow(bpskObj.matchedFiltering(real(nt)));
            
            fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));
            
            csnr = csnrH;
            
            target_csnr(:,1) = csnr_vector;
            for pos = 1 : length(target_csnr)
                
                target = target_csnr(pos, 1);
                if (target == 0)
                    target = 1E-16;
                end
                
                target_csnr(pos, 2) = k*generateK(target, csnr);
            end
            
        end
        
        function [targetCSNR, alpha, khist, performance] = runBPSK()
            %RUNBPSK Application for process BPSK received raw data from the USRP
            %
            
            % getKVector return the vector of k values for the given CSNR
            % assuming a model y = alpha*k*x + n
            %
            
            % output
            
            % params
%             path = '~/Dropbox/Research/data/rx/12-28-12-3/';
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            ks = ones(1, 101); %importdata(strcat(self.sigpath, 'ks/ks101.mat'));
            a_k = importdata(strcat(self.sigpath, self.data));
            
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
       
        function [performance] = runBPSK2(self, isLargeCarrier, show)
            % getKVector return the vector of k values for the given CSNR
            % assuming a model y = alpha*k*x + n
            %
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
            a_k = importdata(strcat(self.sigpath, self.data));

            % instance variables
            bpskObj = ModUtils(ps, self.N, []);
            header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
            
            bpskObj = ModUtils(ps, self.N, ks);
            
            % results
            performance = zeros(1, 3); % Es/N0 BER BERThe
            
            % processing
            pos = 1;
            files = dir(self.pathrx);
            for i = length(files):-1:1
                
                if files(i).isdir == 0 && isempty(strfind(files(i).name, 'raw')) &&...
                        isempty(strfind(files(i).name, 'performance')) &&...
                        isempty(strfind(files(i).name, '.DS_Store'))
                    
                    file = files(i).name;
                    fprintf('i: file = %s\n', strcat(self.pathrx, file));
                    
                    % parameters related to model rt = k*fading*xt + nt
                    
                    % k
                    init = strfind(file, 'at');
                    final = strfind(file, '_G');
                    k = str2double(file(init+2:final-1));
                    
                    % signal name
                    signame = file(1:init-1);
                    dt = read_complex_binary(strcat(self.pathrx, file));
                    
                    % sine
                    sine = read_complex_binary(strcat(self.pathrx, 'fading/', signame, 'at1.0_G0.dat'));
                    
                    % expected nframes to receive
                    lframe = 10000;
                    exp_nframes = floor(length(a_k)/lframe);
                    rec_frames = 0;
                    
                    if (lframe*exp_nframes < length(a_k))
                        exp_nframes = exp_nframes + 1;
                    end
                    
%                     % instance variables
%                     bpskObj = ModUtils(ps, self.N, ks);
%                     header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
                    
                    if isLargeCarrier == true
                        header = header + 1;
                    end
                    smRec = SmartReceiver(header);
                    
                    % Fading estimate fading, assume fading cyclostationary
                    fading = fadingEstimator(sine);
                    
                    fprintf('i: fading = %f, k = %f\n', fading, k);
                    
                    dt = dt(10000:end);
                    
%                     figure; plot(abs(dt));
                    
                    nframes = floor( length(dt) / ((lframe+length(ks) + self.n_preamble + 10)*self.N) );
                    if (((lframe + length(ks) + self.n_preamble + 10)*self.N ) * nframes < length(dt))
                        nframes = nframes+1;
                    end
                    
                    isFirst = true;
                    fixPos = 0;
                    
                    recv = [];
                    nt = [];
                    
                    for fr = 1:nframes
                        
                        % receiving the portion of the message
                        len = (lframe + length(ks) + self.n_preamble + self.gd*4)*self.N;
                        
                        if exp_nframes == 1
                            len2 = ((length(a_k) - lframe*rec_frames) + length(ks) + self.n_preamble + self.gd*4)*self.N;
                            rt = dt(1 + (fr-1)*len + fixPos : 1 + (fr-1)*len + fixPos + len2);
                        elseif exp_nframes > 1
                            rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                        else
                            break;
                        end
                        
                        % rx
                        
                        % timing sync
                        % finding the complete pilot, timing sync
                        if exp_nframes == 1
                            [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), (length(a_k) - lframe*rec_frames)*self.N);
                        else
                            [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), lframe*self.N);
                        end
                        
                        if isempty(headstart) || length(headstart) > 10
                            continue;
                        else
                            
                            wasFound = false;
                            for j = length(headstart):-1:1
                                hs = headstart(j);
                                
                                if isFirst == true
                                    
                                    fprintf('i: is first\n');
                                    isFirst = false;
                                    fixPos = hs - self.n_preamble*self.N;
                                    %                             fixPos = hs - ((self.gd + 1) + self.n_preamble)*self.N;
                                    rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos + 1);
                                    
                                    hs = self.n_preamble*self.N + 1;
                                    
                                end
                                
                                if hs - self.n_preamble*self.N == -1
                                    hs = hs + 1;
                                end
                                
                                if hs >= self.n_preamble*self.N
                                    fprintf('i: headstart found at sample %i', hs);
                                    headstart = hs;
                                    loc = loc(j);
                                    wasFound = true;
                                    
                                    exp_nframes = exp_nframes - 1;
                                    rec_frames = rec_frames + 1;
                                    break;
                                    
                                end
                                
                            end
                            
                            if wasFound == false || length(rt) < headstart+(length(lframe)+length(ks))*self.N;
                                disp('i: a_k not used due to ')
                                if isFirst == false
                                    exp_nframes = exp_nframes - 1;
                                    rec_frames = rec_frames + 1;
                                end
                                continue;
                            end
                            
                            if show == true
                                figure('Name', 'r(t)');
                                subplot(2, 1, 1); hold on; plot(abs(rt)); stem(headstart, abs(rt(headstart)), 'r'); title('r(t)'); hold off;
                                subplot(2, 1, 2); hold on; plot(ycorr); stem(loc, ycorr(loc), 'g'); title('xcorr')
                            end
                            
                        end
                        
                        disp('i: finding marker sequences - OK')
                        
                        %Remove unwanted portions(first few samples till the peak
                        %value)
                        % setting starting sampling point
                        %                 noise = rt(headstart - self.n_preamble*self.N + 1:headstart);
                        %                 rt = rt(init:end);
                        
                        initRt = headstart;
                        initNoise = initRt - self.n_preamble*self.N;
                        
                        noise = rt(initNoise:initNoise + (self.n_preamble - 2*self.gd)*self.N);
                        
                        % larger carrier
                        %                 if isLargeCarrier == true
                        %                     initRt = initRt - (self.gd - 1)*self.N;
                        %                 end
%                         lenPS = length(-self.gd:1/self.N:self.gd);
%                         initRt = initRt - lenPS;
                        
                        % fixing phase
                        if ycorr(loc) < 0
                            disp('i: fixing phase by multiplaying by cos(pi)');
                            rt = rt*exp(1j*pi);
                        end
                        
                        % Getting a more precise init of the Marker Sequence
%                         [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(...
%                             abs(rt(initRt:initRt + 4*length(ks)*self.N)), 0);
                        
%                         initRt = initRt - 1 - (self.gd-2)*self.N;
%                         initRt = initRt + headstart - 1 - (self.gd-2)*self.N;
                        %                 rt = rt(initRt:end);
                        
                        % improving the performance of the frame
                        if length(a_k) > lframe
                            [rk, initRt] = self.improveOutput(bpskObj, initRt, rt, ...
                                [ks a_k((fr-1)*lframe + 1 : fr*lframe)]', k*fading, isLargeCarrier);
                        else
                            [rk, initRt] = self.improveOutput(bpskObj, initRt, rt, ...
                                [ks a_k]', k*fading, isLargeCarrier);
                        end
                        
                        
                        %                 %% fixing phase
                        %                 if ycorr(loc) < 0
                        %                     disp('i: fixing phase by multiplaying by cos(pi)');
                        %                     rt = rt*exp(1j*pi);
                        %                 end
                        %
                        %                 % dem
                        %                 % output of the matched filter and sampling
                        %
                        %                 % removing large carrier
                        %                 if isLargeCarrier == true
                        %                     rk = bpskObj.matchedFiltering(abs(rt) - k*fading);
                        %                 else
                        %                     rk = bpskObj.matchedFiltering(abs(rt));
                        %                 end
                        
                        % header
                        rec_header = bpskObj.decision(rk(1:length(ks)));
                        
                        % message
                        if exp_nframes == 0
                            rec_sig = rk(length(ks) + 1 : length(ks) + 1 + (length(a_k)-lframe*(rec_frames-1))-1);
                        else
                            rec_sig = rk(length(ks) + 1 : length(ks) + 1 + lframe-1);
                        end
                        
                        recv = [recv rec_sig'];
                        nt = [nt noise'];
                        
                    end
                    
                    %% SDR Real test
                    fprintf('\ni: SDR - Testing the output MF and fading*x + n\n');
                    
                    % Ss the theory said
                    
                    % removing large carrier
                    if isLargeCarrier == true
                        vn = var(abs(nt) - k*fading);
                        mn = mean(abs(nt) - k*fading);
                    else
                        vn = var(real(nt));
                        mn = mean(real(nt));
                    end
                    
                    fprintf('i: vn = %d, mn = %d\n', vn, mn);
                    
                    % Performance
                    % CSNR
                    csnr = (pow(recv) / pow(bpskObj.matchedFiltering(real(nt)))) - 1;
                    csnrH = fading^2*k^2*pow(a_k)/pow(bpskObj.matchedFiltering(real(nt)));
                    
                    fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));
                    
                    csnr = csnrH;
                    
                    result = bpskObj.decision(recv);
                    %
                    [number_of_errors, bit_error_rate, pos_errors] = biterr(a_k, result);
                    theoryBer = (1/2)*erfc(sqrt(2*csnr)/sqrt(2));
                    fprintf('i: theory = %d, sim = %d, numerror = %i\n', theoryBer, bit_error_rate, number_of_errors);
                    
                    % Es/N0 BER BERThe
                    performance(pos, 1) = 10*log10(csnr);
                    performance(pos, 2) = bit_error_rate;
                    performance(pos, 3) = theoryBer;
                    pos = pos + 1;
                    
                    
                end
            end
            
            % Plotting
            performance = sortrows(performance, -1);
            
            simBer = performance(:, 2); % simulated ber
            theoryBer = performance(:, 3); % theoretical ber
            
            figure;
            semilogy(performance(:, 1), theoryBer','b.-');
            hold on;
            semilogy(performance(:, 1), simBer','mx-');
            hold off;
            
            grid on
            legend('theory', 'real transmission');
            xlabel('Eb/No, dB');
            ylabel('Bit Error Rate');
            title('testBERPilotPlots - Bit error probability curve for BPSK modulation');
            
        end
        
        function [performance] = runBPSK3(self, isLargeCarrier, show)
            % getKVector return the vector of k values for the given CSNR
            % assuming a model y = alpha*k*x + n
            %
            
            % signals
            ps = importdata(strcat(self.sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', self.N, self.gd, self.ro))); 
            ps = ps/sqrt(sum((ps).^2));
            
            ks = importdata(strcat(self.sigpath, 'ks/', self.ksf));
            a_k = importdata(strcat(self.sigpath, self.data));

            % instance variables
            bpskObj = ModUtils(ps, self.N, []);
            header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
            
            bpskObj = ModUtils(ps, self.N, ks);
            
            % results
            performance = zeros(1, 3); % Es/N0 BER BERThe
            
            % processing
            pos = 1;
            files = dir(self.pathrx);
            for i = length(files):-1:1
                
                if files(i).isdir == 0 && isempty(strfind(files(i).name, 'raw')) &&...
                        isempty(strfind(files(i).name, 'performance')) &&...
                        isempty(strfind(files(i).name, '.DS_Store'))
                    
                    file = files(i).name;
                    fprintf('i: file = %s\n', strcat(self.pathrx, file));
                    
                    % parameters related to model rt = k*fading*xt + nt
                    
                    % k
                    init = strfind(file, 'at');
                    final = strfind(file, '_G');
                    k = str2double(file(init+2:final-1));
                    
                    % signal name
                    signame = file(1:init-1);
                    dt = read_complex_binary(strcat(self.pathrx, file));
                    
                    % sine
                    sine = read_complex_binary(strcat(self.pathrx, 'fading/', signame, 'at1.0_G0.dat'));
                    
                    % expected nframes to receive
                    lframe = 10000;
                    exp_nframes = floor(length(a_k)/lframe);
                    rec_frames = 0;
                    
                    if (lframe*exp_nframes < length(a_k))
                        exp_nframes = exp_nframes + 1;
                    end
                    
%                     % instance variables
%                     bpskObj = ModUtils(ps, self.N, ks);
%                     header =  bpskObj.pulseShaping(bpskObj.bpskMod(ks, 0, 0)');
                    
                    if isLargeCarrier == true
                        header = header + 1;
                    end
                    smRec = SmartReceiver(header);
                    
                    % Fading estimate fading, assume fading cyclostationary
                    fading = fadingEstimator(sine);
                    
                    fprintf('i: fading = %f, k = %f\n', fading, k);
                    
                    dt = dt(10000:end);
                    
                    nframes = floor( length(dt) / ((lframe+length(ks) + self.n_preamble + 10)*self.N) );
                    if (((lframe + length(ks) + self.n_preamble + 10)*self.N ) * nframes < length(dt))
                        nframes = nframes+1;
                    end
                    
                    isFirst = true;
                    fixPos = 0;
                    
                    recv = [];
                    nt = [];
                    
                    for fr = 1:nframes
                        
                        % receiving the portion of the message
                        len = (lframe + length(ks) + self.n_preamble + self.gd*4)*self.N;
                        
                        if exp_nframes == 1
                            len2 = ((length(a_k) - lframe*rec_frames) + length(ks) + self.n_preamble + self.gd*4)*self.N;
                            rt = dt(1 + (fr-1)*len + fixPos : 1 + (fr-1)*len + fixPos + len2);
                        elseif exp_nframes > 1
                            rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos);
                        else
                            break;
                        end
                        
                        % rx
                        
                        % timing sync
                        % finding the complete pilot, timing sync
                        if exp_nframes == 1
                            [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), (length(a_k) - lframe*rec_frames)*self.N);
                        else
                            [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(abs(rt), lframe*self.N);
                        end
                        
                        if isempty(headstart) || length(headstart) > 10
                            continue;
                        else
                            
                            wasFound = false;
                            for j = length(headstart):-1:1
                                hs = headstart(j);
                                
                                if isFirst == true
                                    
                                    fprintf('i: is first\n');
                                    isFirst = false;
                                    fixPos = hs - self.n_preamble*self.N;
                                    %                             fixPos = hs - ((self.gd + 1) + self.n_preamble)*self.N;
                                    rt = dt(1 + (fr-1)*len + fixPos : fr*len + fixPos + 1);
                                    
                                    hs = self.n_preamble*self.N + 1;
                                    
                                end
                                
                                if hs - self.n_preamble*self.N == -1
                                    hs = hs + 1;
                                end
                                
                                if hs >= self.n_preamble*self.N
                                    fprintf('i: headstart found at sample %i', hs);
                                    headstart = hs;
                                    loc = loc(j);
                                    wasFound = true;
                                    
                                    exp_nframes = exp_nframes - 1;
                                    rec_frames = rec_frames + 1;
                                    break;
                                    
                                end
                                
                            end
                            
                            if wasFound == false || length(rt) < headstart+(length(lframe)+length(ks))*self.N;
                                disp('i: a_k not used due to ')
                                if isFirst == false
                                    exp_nframes = exp_nframes - 1;
                                    rec_frames = rec_frames + 1;
                                end
                                continue;
                            end
                            
                            if show == true
                                figure('Name', 'r(t)');
                                subplot(2, 1, 1); hold on; plot(abs(rt)); stem(headstart, abs(rt(headstart)), 'r'); title('r(t)'); hold off;
                                subplot(2, 1, 2); hold on; plot(ycorr); stem(loc, ycorr(loc), 'g'); title('xcorr')
                            end
                            
                        end
                        
                        disp('i: finding marker sequences - OK')
                        
                        %Remove unwanted portions(first few samples till the peak
                        %value)
                        % setting starting sampling point
                        %                 noise = rt(headstart - self.n_preamble*self.N + 1:headstart);
                        %                 rt = rt(init:end);
                        
                        initRt = headstart;
                        initNoise = initRt - self.n_preamble*self.N;
                        
                        noise = rt(initNoise:initNoise + (self.n_preamble - 2*self.gd)*self.N);
                        
                        % larger carrier
                        %                 if isLargeCarrier == true
                        %                     initRt = initRt - (self.gd - 1)*self.N;
                        %                 end
                        lenPS = length(-self.gd:1/self.N:self.gd);
                        initRt = initRt - lenPS;
                        
                        % fixing phase
                        if ycorr(loc) < 0
                            disp('i: fixing phase by multiplaying by cos(pi)');
                            rt = rt*exp(1j*pi);
                        end
                        
                        % Getting a more precise init of the Marker Sequence
                        [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence2(...
                            abs(rt(initRt:initRt + 4*length(ks)*self.N)), 0);
                        
                        initRt = initRt + headstart - 1 - (self.gd-2)*self.N;
                        %                 rt = rt(initRt:end);
                        
                        % improving the performance of the frame
                        if length(a_k) > lframe
                            [rk, initRt] = self.improveOutput(bpskObj, initRt, rt, ...
                                [ks a_k((fr-1)*lframe + 1 : fr*lframe)]', k*fading, isLargeCarrier);
                        else
                            [rk, initRt] = self.improveOutput(bpskObj, initRt, rt, ...
                                [ks a_k]', k*fading, isLargeCarrier);
                        end
                        
                        
                        %                 %% fixing phase
                        %                 if ycorr(loc) < 0
                        %                     disp('i: fixing phase by multiplaying by cos(pi)');
                        %                     rt = rt*exp(1j*pi);
                        %                 end
                        %
                        %                 % dem
                        %                 % output of the matched filter and sampling
                        %
                        %                 % removing large carrier
                        %                 if isLargeCarrier == true
                        %                     rk = bpskObj.matchedFiltering(abs(rt) - k*fading);
                        %                 else
                        %                     rk = bpskObj.matchedFiltering(abs(rt));
                        %                 end
                        
                        % header
                        rec_header = bpskObj.decision(rk(1:length(ks)));
                        
                        % message
                        if exp_nframes == 0
                            rec_sig = rk(length(ks) + 1 : length(ks) + 1 + (length(a_k)-lframe*(rec_frames-1))-1);
                        else
                            rec_sig = rk(length(ks) + 1 : length(ks) + 1 + lframe-1);
                        end
                        
                        recv = [recv rec_sig'];
                        nt = [nt noise'];
                        
                    end
                    
                    %% SDR Real test
                    fprintf('\ni: SDR - Testing the output MF and fading*x + n\n');
                    
                    % Ss the theory said
                    
                    % removing large carrier
                    if isLargeCarrier == true
                        vn = var(abs(nt) - k*fading);
                        mn = mean(abs(nt) - k*fading);
                    else
                        vn = var(real(nt));
                        mn = mean(real(nt));
                    end
                    
                    fprintf('i: vn = %d, mn = %d\n', vn, mn);
                    
                    % Performance
                    % CSNR
                    csnr = (pow(recv) / pow(bpskObj.matchedFiltering(real(nt)))) - 1;
                    csnrH = fading^2*k^2*pow(a_k)/pow(bpskObj.matchedFiltering(real(nt)));
                    
                    fprintf('i: csnr: %12.5f; csnrH: %12.5f\n', 10*log10(csnr), 10*log10(csnrH));
                    
                    csnr = csnrH;
                    
                    result = bpskObj.decision(recv);
                    %
                    [number_of_errors, bit_error_rate, pos_errors] = biterr(a_k, result);
                    theoryBer = (1/2)*erfc(sqrt(2*csnr)/sqrt(2));
                    fprintf('i: theory = %d, sim = %d, numerror = %i\n', theoryBer, bit_error_rate, number_of_errors);
                    
                    % Es/N0 BER BERThe
                    performance(pos, 1) = 10*log10(csnr);
                    performance(pos, 2) = bit_error_rate;
                    performance(pos, 3) = theoryBer;
                    pos = pos + 1;
                    
                    
                end
            end
            
            % Plotting
            performance = sortrows(performance, -1);
            
            simBer = performance(:, 2); % simulated ber
            theoryBer = performance(:, 3); % theoretical ber
            
            figure;
            semilogy(performance(:, 1), theoryBer','b.-');
            hold on;
            semilogy(performance(:, 1), simBer','mx-');
            hold off;
            
            grid on
            legend('theory', 'real transmission');
            xlabel('Eb/No, dB');
            ylabel('Bit Error Rate');
            title('testBERPilotPlots - Bit error probability curve for BPSK modulation');
            
        end
    
        
    end
    
    %% Static methods
    methods(Static)
        function main(option, show)
            % MAIN run the desired procedure.
            % MAIN(option, args). The options are:
            %   1: to run the creation of a BPSK signal
            %
            
            %             if nargin < 5, subdiv = 20; end
            %             if nargin < 4, angle = 10; end
            %             if nargin < 3, npts = 25; end
            if nargin < 2, show = false; end
            
            
            switch option
                case 1
                    disp('i: create BPSK signal')
                    
                    N = 20;
                    gd = 5;
                    ro = 20;
                    sigpath = 'data/';
                    data = 'a_k/a_k1E4.mat';
                    np = 500;
                    pathto = '~/Dropbox/Research/data/tx/default/';
                    ksf = 'ks500ones.mat';
                    lc = false;
                    
                    app = AppBPSK(N, gd, ro, sigpath, data, ksf, np, '');

                    st = app.createBSPKSignal(lc, 'bpsk1E4ks500big.dat', pathto);
                    
                    if show == true
                       figure; plot(st);
                    end
                    
                    
                case 2
                    disp('i: get K attenuation factor for transmission')
                    
                    N = 20;
                    gd = 5;
                    ro = 20;
                    sigpath = 'data/';
                    data = 'a_k/a_k1E4.mat';
                    np = 500;
                    pathrx = '~/Dropbox/Research/data/rx/02-05-13/';
                    ksf = 'ks500ones.mat';
                    lc = true;
                    
                    app = AppBPSK(N, gd, ro, sigpath, data, ksf, np, pathrx);

                    target_csnr = app.getKVector2(0:5:30, 'bpsklc1E4ks500v2_at0.5_G0.dat', lc, show);
                    
                    disp(sprintf('%f;', target_csnr(:, 2)));
                    
                    
                case 3
                    disp('i: run BPSK')
                    
                    N = 20;
                    gd = 5;
                    ro = 20;
                    sigpath = 'data/';
                    data = 'a_k/a_k1E4.mat';
                    np = 500;
                    pathrx = '~/Dropbox/Research/data/rx/02-05-13/';
                    ksf = 'ks500ones.mat';
                    lc = true;
                    
                    app = AppBPSK(N, gd, ro, sigpath, data, ksf, np, pathrx);

                    performance = app.runBPSK2(lc, show)
                    
                    
                case 4
                    disp('i: run BPSK3')
                    
                    N = 20;
                    gd = 5;
                    ro = 20;
                    sigpath = 'data/';
                    data = 'a_k/a_k1E4.mat';
                    np = 500;
                    pathrx = '~/Dropbox/Research/data/rx/01-16-13-1/';
                    ksf = 'ks500ones.mat';
                    
                    app = AppBPSK(N, gd, ro, sigpath, data, ksf, np, pathrx);
                    
                    performance = app.runBPSK3(true, show)
                    
                case 5
                    disp('i: create 2 BPSK signal')
                    
                    N = 20;
                    gd = 5;
                    ro = 20;
                    sigpath = 'data/';
                    data = 'a_k/a_k1E4.mat';
                    np = 500;
                    pathto = '~/Dropbox/Research/data/tx/default/';
                    ksf = 'ks500ones.mat';
                    lc = true;
                    
                    app = AppBPSK(N, gd, ro, sigpath, data, ksf, np, '');
                    
                    st = app.createBSPKSignal2('bpsklc1E4ks500v2.dat', pathto);
                    
                    if show == true
                        figure; plot(st);
                    end
                    
                otherwise
                    disp('Unknown method.')
            end
            
        end
    end
end

