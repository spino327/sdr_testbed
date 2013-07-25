%Copyright (c) 2011, Universidad Industrial de Santander, Colombia
%University of Delaware
%All rights reserved.
%@author: Sergio Pino
%@author: Henry Arguello
%Website: http://www.eecis.udel.edu/
%emails  : sergiop@udel.edu - henarfu@udel.edu
%Date   : March, 2011

classdef SmartReceiver
    %SMARTRECEIVER This is the Smart receiver implementation
    
    properties
        marker_seq
    end
    
    methods
        
        function obj = SmartReceiver(ms)
        % constructor obj = SmartReceiver(ms)
        % in:
        %   - ms = known sequence or market sequence
        % out:
        %   - SmartReceiver object
            if size(ms, 1) > 1
                ms = ms';
            end
            
            obj.marker_seq = ms; 
        end
        
        function [headstart, corr, peaks, loc] = findMarkerSequence(obj, data, len)
        % findMarkerSequence [headstart, ycorr, peaks, loc] = findMarkerSequence(obj, data)
        % in:
        %   - data = samples data
        %   - len = length message
        % out:
        %   - headstart = sampling start time for each Marker Sequence find
        %   in data
        %   - ycorr = correlation result sequence, this is between data and
        %   marker sequence in the ycorr sequence
        %   - peaks = correlation peak value
        %   - loc = peak localization in the ycorr sequence
        
        
        % fixme sergiop@udel.edu we need to work recursively
        
            % find problematic peaks 
            fixData = abs(data);
            meanData = mean(fixData);
            varData = var(fixData);
            [peaks, loc] = findpeaks(fixData, 'MINPEAKHEIGHT', 3*sqrt(varData) + 3*meanData);
            
            for i=1:length(loc)
                
                % killing left part
                j = 1;
                exit = false;
                while loc(i)-j >= 1 && exit ~= true
                    if fixData(loc(i)-j) > 3*sqrt(varData) + 3*meanData
                        data(loc(i)-j) = 0;
                    else
                        exit = true;
                    end
                    j = j + 1;
                end
                
                % killing right part
                j = 1;
                exit = false;
                while loc(i)+j <= length(data) && exit ~= true
                    if fixData(loc(i)+j) > 3*sqrt(varData) + 3*meanData
                        data(loc(i)+j) = 0;
                    else
                        exit = true;
                    end
                    j = j + 1;
                end
                
            end
            
            if size(data, 1) > 1
                data = data';
            end
            
%             corr = fliplr(ifft(fft(fliplr(data)).*fft([obj.marker_seq zeros(1, length(data)-length(obj.marker_seq))])));
%             ycorr = abs(corr);
            
            corr = xcorr(obj.marker_seq, data);
            ycorr = abs(corr);           % do cross correlation
            
            k = sqrt(ycorr*ycorr'/length(ycorr));
            ycorr = ycorr/k;
            
            [peaks, loc] = findpeaks(ycorr, 'MINPEAKHEIGHT', 0.80*max(ycorr), 'MINPEAKDISTANCE', len);%length(obj.marker_seq) + len);
            headstart = length(data)-loc;            % gives place where header starts
            
            for h=1:length(headstart)
                
                if headstart(h) <= 0
                    
                    if h == 1
                        headstart = headstart(h+1:end);
                    
                    elseif h == length(headstart)
                        headstart = headstart(1:end-1);
                    
                    else
                        headstart = [headstart(1:h-1) headstart(h+1:end)];
                        
                    end
                end
                
            end
            
        end
        
        function [headstart, corr, peaks, loc] = findMarkerSequence2(obj, data, len)
        % findMarkerSequence [headstart, ycorr, peaks, loc] = findMarkerSequence(obj, data)
        % in:
        %   - data = samples data
        %   - len = length message
        % out:
        %   - headstart = sampling start time for each Marker Sequence find
        %   in data
        %   - ycorr = correlation result sequence, this is between data and
        %   marker sequence in the ycorr sequence
        %   - peaks = correlation peak value
        %   - loc = peak localization in the ycorr sequence
        
        
        % fixme sergiop@udel.edu we need to work recursively
        
            % find problematic peaks 
%             fixData = abs(data);
%             meanData = mean(fixData);
%             varData = var(fixData);
%             [peaks, loc] = findpeaks(fixData, 'MINPEAKHEIGHT', 3*sqrt(varData) + 3*meanData);
%             
%             for i=1:length(loc)
%                 
%                 % killing left part
%                 j = 1;
%                 exit = false;
%                 while loc(i)-j >= 1 && exit ~= true
%                     if fixData(loc(i)-j) > 3*sqrt(varData) + 3*meanData
%                         data(loc(i)-j) = 0;
%                     else
%                         exit = true;
%                     end
%                     j = j + 1;
%                 end
%                 
%                 % killing right part
%                 j = 1;
%                 exit = false;
%                 while loc(i)+j <= length(data) && exit ~= true
%                     if fixData(loc(i)+j) > 3*sqrt(varData) + 3*meanData
%                         data(loc(i)+j) = 0;
%                     else
%                         exit = true;
%                     end
%                     j = j + 1;
%                 end
%                 
%             end
            
            if size(data, 1) > 1
                data = data';
            end
            
            corr = xcorr(obj.marker_seq, data);
            ycorr = abs(corr);           % do cross correlation
            
            k = sqrt(ycorr*ycorr'/length(ycorr));
            ycorr = ycorr/k;
            
            [peaks, loc] = max(ycorr);

            headstart = length(data)-loc;            % gives place where header starts
            
            for h=1:length(headstart)
                
                if headstart(h) <= 0
                    
                    if h == 1
                        headstart = headstart(h+1:end);
                    
                    elseif h == length(headstart)
                        headstart = headstart(1:end-1);
                    
                    else
                        headstart = [headstart(1:h-1) headstart(h+1:end)];
                        
                    end
                end
                
            end
            
        end
        
        function [headstart, corr, peaks, loc] = findSequence(obj, data, seq, mpd, np)
        % FINDSEQUENCE [headstart, ycorr, peaks, loc] = FINDSEQUENCE(data, mpd, np)
        % in:
        %   - data = samples data
        %   - seq = sequence to find
        %   - mpd = length message, 'minpeakdistance' returns only peaks
        %   with indices separated by more than the positive integer mpd.
        %   mpd defaults to 1.
        %   - np = returns a maximum number of peaks np.
        % out:
        %   - headstart = sampling start time for each Marker Sequence find
        %   in data
        %   - ycorr = correlation result sequence, this is between data and
        %   marker sequence in the ycorr sequence
        %   - peaks = correlation peak value
        %   - loc = peak localization in the ycorr sequence
        
        
        % fixme sergiop@udel.edu we need to work recursively
        
            % find problematic peaks 
            fixData = abs(data);
            meanData = mean(fixData);
            varData = var(fixData);
            [peaks, loc] = findpeaks(fixData, 'MINPEAKHEIGHT', 3*sqrt(varData) + 3*meanData);
            
            for i=1:length(loc)
                
                % killing left part
                j = 1;
                exit = false;
                while loc(i)-j >= 1 && exit ~= true
                    if fixData(loc(i)-j) > 3*sqrt(varData) + 3*meanData
                        data(loc(i)-j) = 0;
                    else
                        exit = true;
                    end
                    j = j + 1;
                end
                
                % killing right part
                j = 1;
                exit = false;
                while loc(i)+j <= length(data) && exit ~= true
                    if fixData(loc(i)+j) > 3*sqrt(varData) + 3*meanData
                        data(loc(i)+j) = 0;
                    else
                        exit = true;
                    end
                    j = j + 1;
                end
                
            end
            
            if size(data, 1) > 1
                data = data';
            end
            
            corr = xcorr(seq, data);
            ycorr = abs(corr);
            k = sqrt(ycorr*ycorr'/length(ycorr));
            ycorr = ycorr/k;
            
            [peaks, loc] = findpeaks(ycorr, 'MINPEAKHEIGHT', 0.80*max(ycorr), 'MINPEAKDISTANCE', mpd);%length(obj.marker_seq) + len);

            headstart = length(data)-loc;            % gives place where header starts
            
            if np ~= 1
            
                for h=1:length(headstart)

                    if headstart(h) <= 0

                        if h == 1
                            headstart = headstart(h+1:end);

                        elseif h == length(headstart)
                            headstart = headstart(1:end-1);

                        else
                            headstart = [headstart(1:h-1) headstart(h+1:end)];

                        end
                    end

                end
                
            else
                
                headstart = headstart(end);
                
            end
            
        end
        
        function df = lR(self, x, N)
            L = length(x);
            R = zeros(1, N);
            for m=1:N
                R(m)=sum(x(m+1:L).*conj(x(1:L-m)))/(L-m);
            end
            m = 1:N;
            
            df = sum(m.*(abs(R).^2).*angle(R)) / sum((m.^2).*(abs(R).^2));
        end
        
        function [df] = fitz (self, x, N)
        % frequency synchronization with fitz algorithm 
            L = length(x);
            R=zeros(1,N);
            for m=1:N
                R(m)=sum(x(m+1:L).*conj(x(1:L-m)))/(L-m);
            end
            %f=abs(sum(angle(R))/pi/N/(N+1)*fmuestreo);
            df=sum(angle(R))/N/(N+1)*2;
        end
        
        function df = kay(self, x)
            L=length(x);
            gamma=1.5.*L./(L.^2-1).*(1-(2.*(1:L-1)./L-1).^2);
            xdesp=x(2:L);
            %     f=abs(sum(gamma.*angle(xdesp.*conj(x(1:L-1))))/2/pi*fmuestreo);
            df = 1/(2*pi)*sum(gamma.*angle(xdesp.*conj(x(1:L-1))));
        end
        
        function [rkFix, cFix] = frequencySync(self, rk)
            
            % setup
            rkFix = rk;
            snr = pow(real(rkFix))/pow(imag(rkFix));
            dt = -1;
 
            status = true;
            i = 0;
            cFix = 0;
            while (status)
                
                i = i + 1;
                
                fix = self.kay(imag(rkFix(1:floor(end/2)))');
                cFix = cFix + fix;
                rkT = rkFix*exp(dt*1i*fix);
%                 rkT = rkFix*exp(dt*1i*self.kay(imag(rkFix(1:floor(end/2)))'));
%                 rkT = rkFix*exp(dt*1i*self.fitz(imag(rkFix)', floor(length(rk)/4)));
%                 rkT = rkFix*exp(dt*1i*self.lR(imag(rkFix(1:floor(end/2)))', floor(length(rk)/4)));
                snrT = pow(real(rkT))/pow(imag(rkT));

                % SNR
                if (snrT > snr)
                    snr = snrT;
                    rkFix = rkT;
                    
                elseif (snrT <= snr) && (i <= 2)
                    
                    dt = dt*-1;
                    cFix = 0;
                    disp('i: freq sync - dir change');

                elseif (i > 2)
%                     disp(sprintf('i: freq sync - in %i', i));
                    break;
                end
                        
            end           
            
            cFix = dt*cFix;
        end
        
        function [recv, nt] = recvPackage(self, dtF, lframe, lxt, exp_nframes, bpskObj, preamble)
            % Receive the packages
            % [recv, nt] = recvPackage(self, dtF, lframe, lxt, exp_nframes, bpskObj, preamble)
            % input:
            %   - dtF = signal in the antenna
            %   - lframe
            %   - lxt = length xt
            %   - exp_nframes
            %   - bpskObj
            %   - preamble
            %
            % output:
            %   - recv = received signal
            %   - nt = noise
            
            rec_frames = 0;
            
            ks = bpskObj.know_seq;
            N = bpskObj.N;
            
            nframes = floor(length(dtF)/((lframe+length(ks)+preamble + 10)*N));
            if (((lframe+length(ks)+preamble + 10)*N)*nframes < length(dtF))
                nframes = nframes+1;
            end
            
            isFirst = true;
            fixPos = 0;
            
            recv = [];
            nt = [];
            
            for fr=1:nframes
                
                % parameters
                diff = 0; % measures the difference between hs - preamble*N;
                
                % receiving the portion of the message
                len = (lframe+length(ks)+preamble + 10)*N;
                rt = dtF(1 + (fr-1)*len + fixPos:fr*len + fixPos);
                
                % matched filtering
                rtf = bpskObj.mFilter(rt);
                [rtf, cfix] = self.frequencySync(rtf);
                rt = rt*exp(1i*cfix);
                fprintf('i: cumulative frequency fix = %f\n', cfix);
                
                %% rx
                
                % timing sync
                % finding the complete pilot, timing sync
                if exp_nframes == 1
                    [headstart, ycorr, peaks, loc] = self.findMarkerSequence(real(rt), (lxt-lframe*rec_frames)*N);
                else
                    [headstart, ycorr, peaks, loc] = self.findMarkerSequence(real(rt), lframe*N);
                end
                
                if isempty(headstart) || length(headstart) > 10
                    continue;
                else
                    
                    wasFound = false;
                    for j = length(headstart):-1:1
                        hs = headstart(j);
                        
                        rtf = bpskObj.mFilter(rt);
                        if (pow(real(rtf))/pow(imag(rtf))-1)>=0.1
                            
                            if isFirst
                                
                                fprintf('i: is first\n');
                                isFirst = false;
                                fixPos = hs - preamble*N;
                                rt = dtF(1 + (fr-1)*len + fixPos:fr*len + fixPos);
                                
                                % matched filtering
                                rtf = bpskObj.mFilter(rt);
                                [rtf, cfix] = self.frequencySync(rtf);
                                rt = rt*exp(1i*cfix);
                                fprintf('i: cumulative frequency fix = %f\n', cfix);
                                
                                hs = preamble*N + 1;
                                
                            end
                            
                            % How many samples we tolerate to lost
                            lost = 5;
                            diff = hs - preamble*N;
%                             if sign(diff) == -1 &&...
%                                     abs(diff) <= lost
%                                 hs = hs + abs(diff);
%                             end
                            if (sign(diff) == 1 ||...
                                    abs(diff) > lost)
                                diff = 0;
                            end
                            
                            if hs + abs(diff) >= preamble*N
                                fprintf('i: headstart found at sample %i\n', hs);
                                headstart = hs;
                                loc = loc(j);
                                wasFound = true;
                                
                                exp_nframes = exp_nframes - 1;
                                rec_frames = rec_frames + 1;
                                break;
                                
                            end
                        end
                    end
                    
                    if (wasFound == false ||...
                            length(rt) < headstart+(lframe+length(ks))*N)
                        
                        disp('i: data not used due to ');
                        if isFirst == false
                            exp_nframes = exp_nframes - 1;
                            rec_frames = rec_frames + 1;
                        end
                        continue;
                    end
                    
                end
                
                disp('i: finding marker sequences - OK');
                
                % setting starting sampling point
%                 init = headstart + 1;
                init = headstart;
                
                %Remove unwanted portions(first few samples till the peak value)
                noise = rt(headstart - preamble*N + abs(diff) + 1:headstart);
                rt = rt(init:end);
                
                %% fixing phase
                if ycorr(loc) < 0
                    disp('i: fixing phase by multiplaying by cos(pi)');
                    rt = rt*cos(pi);
                end
                
                % dem
                rk = bpskObj.matchedFiltering(real(rt)); % output of the matched filter and sampling
                
                % header
                result = bpskObj.decision(rk(1:length(ks)));
                
                % message
                if exp_nframes == 0
                    result = rk(length(ks) + 1:length(ks) + 1 + (lxt - lframe*(rec_frames-1))-1);
                else
                    result = rk(length(ks) + 1:length(ks) + 1 + lframe-1);
                end

                recv = [recv result'];
                nt = [nt noise'];
                
                if exp_nframes == 0
                    break;
                end
                
            end
            
        end
          
    end
    
    methods(Static)
        function main(option, show)
            % MAIN run the desired procedure.
            % MAIN(option, args). The options are:
            %   
            
            if nargin < 2, show = false; end
            
            
            switch option
                case 1
                    disp('i: find optimal pilot sequence lenght given the snr')
                    
                    step = 50;
                    NS = 20;
                    gd = 5;
                    ro = 20;
                    sigpath = 'data/';
                    data = 'a_k/a_k1E4.mat';
                    np = 500;
                    pathto = 'data/kx/';
                                        
                    CSNR = 20;
                    sigman=1/sqrt(10^(0.1*CSNR));
                    
                    % signals
                    ps = importdata(strcat(sigpath, sprintf('rcosine/srrc_up%i_gd%i_ro%i.mat', NS, gd, ro)));
                    ps = ps/sqrt(sum((ps).^2));
                    
                    a_k = importdata(strcat(sigpath, data));
                    
                    headerObj = ModUtils(ps, NS, []);
                    
                    % optimizing
                    output = zeros(length(100:50:2000), 2);
                    
                    pos = 1;
                    for len=100:step:2000
                        
                        ks = ones(1, len);
                        % ks
                        for i=1:len/step
                            ks(1 + (i-1)*step : i*step) = mod(i, 2);
                        end
                    
                        % receiver
                        smRec = SmartReceiver(headerObj.pulseShaping(headerObj.bpskMod(ks, 0, 0)));
                        
                        % mod
                        bpskObj = ModUtils(ps, NS, ks);
                        [sym] = bpskObj.bpskMod(a_k, np, 0);
                        st = bpskObj.pulseShaping(sym);
                        
                        % noise
                        noise = sigman*randn(1, (np + len + length(a_k))*NS + length(ps) - 1);
                        
                        % y
                        rt = noise.*st;
                        
                        % signal to noise ratio
                        snr = pow(a_k)/pow(bpskObj.matchedFiltering(noise));
                        
                        % getting the signal
                        headstart = smRec.findMarkerSequence2(rt, 1);
                        
                        recv = bpskObj.matchedFiltering(rt(headstart:end));
                        result = bpskObj.decision(recv);
                        
                        expected = [ks a_k];
                        [number_of_errors, bit_error_rate, pos_errors] = biterr(expected(1:length(result)), result);
%                         theoryBer = (1/2)*erfc(sqrt(2*csnr)/sqrt(2));
                        
                        output(pos, 1) = len;
                        output(pos, 2) = bit_error_rate;
                        
                        pos = pos + 1;
                    end
                                        
                    maxsdr = max(output);

                    
                case 2
                    
                    step = 50;
                    N = 20;
                    np = 1000;
                    
                    NS = 10^6; % number of bits or symbols
                    rand('state',100); % initializing the rand() function
                    randn('state',200); % initializing the randn() function
                    
                    % Transmitter
                    ip = rand(1,NS)>0.5; % generating 0,1 with equal probability
%                     s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 1
                    n = 1/sqrt(2)*[randn(1,NS) + j*randn(1,NS)]; % white gaussian noise, 0dB variance
                    Eb_N0_dB = 10;%0:5:30; % multiple Eb/N0 values
                    
                    headerObj = ModUtils([], NS, []);
                    
                    for ii = 1:length(Eb_N0_dB)
                        
                        for len=100:step:2000
                            
                            ks = ones(1, len);
                            % ks
                            for i=1:len/step
                                ks(1 + (i-1)*step : i*step) = mod(i, 2);
                            end
                            
                            % receiver
                            smRec = SmartReceiver(headerObj.bpskMod(ks, 0, 0));
                            
                            % mod
                            bpskObj = ModUtils([], N, ks);
                            [s] = bpskObj.bpskMod(ip, np, 0);
                            
                            % Noise addition
                            
                            noise = 10^(-Eb_N0_dB(ii)/20)*[n 1/sqrt(2)*[randn(1, len + np) + j*randn(1, len + np)]]; 
                            y = s + noise;
                            
                            % getting the signal
                            headstart = smRec.findMarkerSequence2(y, 1);
                            
                            % receiver - hard decision decoding
                            ipHat = real(y(headstart+len:end))>0;
                            
                            % counting the errors
                            nErr(ii) = size(find([ip - ipHat]),2);
                            
                            % simulated ber
                            simBer = nErr/NS; 
                        end

                    end
                    
                otherwise
                    disp('Unknown method.')
            end
            
        end
    end
    
end

