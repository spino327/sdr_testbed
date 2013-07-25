%Copyright (c) 2011, Universidad Industrial de Santander, Colombia
%University of Delaware
%All rights reserved.
%@author: Sergio Pino
%@author: Henry Arguello
%Website: http://www.eecis.udel.edu/
%emails  : sergiop@udel.edu - henarfu@udel.edu
%Date   : Jun, 2011

classdef ModUtils
    %ModUtils This class represents the modulations and demodulation
    %   - 
    %   -
    %   -
    %   -  Detailed explanation goes here
    
    properties
        %   pulseShape = sampling data of the desired pulse shape
        pulseShape
        %   N = Upsampling factor of the pulse shape (samples per period)
        N
        % know sequence
        know_seq = [1 0];
        % Limiting the response to -GDT to GDT
        % gd = Number of symbol periods between the start of the Raised-cosine and its peak
        GD
    end
    
    methods
        
        function obj = ModUtils(ps, N, ks)
            % constructor for ModUtils objects obj = ModUtils(ps, N, ks)
            %
            % in:   
            %   - ps = pulse shape
            %   - N = upsampling factor
            %   - ks = known sequence PN
            % out: BPSK Object
            if size(ps, 1) > 1
               ps = ps'; 
            end
            
            %Normalization to unit energy
%             ps = ps/sqrt(ps*ps'/length(ps));
            obj.pulseShape = ps;
            
            obj.N = N;
            
            if size(ks, 1) > 1
               ks = ks'; 
            end
            
            obj.know_seq = ks;
            obj.GD = (length(obj.pulseShape) - 1)/(2*obj.N);
        end
        
        function [st, es] = mod(obj, a_k, n_preamble, pos_preamble, debug)
            
            %mod Perform the BPSK modulation using the given parameters: [st, es] = mod(obj, a_k, n_preamble, pos_preamble, debug)
            %   v0.0.1   
            %
            %   [st, es] = mod(obj, a_k, n_preamble, pos_preamble, debug)
            %   
            %   in:
            %
            %   a_k = data vector (only 1 or 0)
            %   n_preamble = number of preamble symbols (0s)
            %   pos_preamble = number of pos preamble symbols (0s)
            %   debug = 1 or 0, if 1 plots at least the first 150 symbols and save a file of st
            %
            %   out:
            %
            %   st = sample modulate signal, [preamble][know seq][five 0s][bpsk mod]
            %   es = pulse shape original enery
            %

%             [es, esd] = signalEnergy(obj.pulseShape, 1, Inf, obj.N, 0);
            es = pow(obj.pulseShape);

            GD = (length(obj.pulseShape) - 1)/(2*obj.N);%  GD = Number of symbol periods between the start of the Raised-cosine and its peak

            st = zeros(length(obj.pulseShape)+obj.N*(length(a_k)-1), 1); % baseband BPSK signal

            disp(sprintf('i:length a_k = %i, length st = %i', length(a_k), length(st)));

            for i=0:length(a_k)-1

                if a_k(i+1) == 1
                    theta_m = 0;
                else     
                    theta_m = pi;
                end

                aigt_iT = obj.pulseShape*cos(theta_m);

                for j=1:length(aigt_iT)

                   st(i*obj.N+j)=st(i*obj.N+j)+aigt_iT(j); 

                end

            end

            % preamble
            if n_preamble > 0
                
                preamble = (1/(max(obj.pulseShape)))*obj.mod(obj.know_seq, 0, 0, 0);

                % preamble of zeros
                preamble = [zeros(obj.N*n_preamble, 1); preamble; zeros(obj.N*pos_preamble, 1)];

                % signal s(t)         
                st = [preamble; st];
                disp(sprintf('i:lenght of st without preamble = %i, st with preamble = %i', length(st)-obj.N*n_preamble, length(st)));

            end

            % debug mode
            if debug == 1

                if length(a_k) > 150
                    toEnd = 150;
                else
                    toEnd = length(a_k);
                end

                % a_k plot 
                figure;
                stem(a_k(1:toEnd));
                title('a_k');

                % s(t) plot 
                figure;
                stem(st(1:(2*obj.N*GD+1)+obj.N*(toEnd-1)), 'rs');
                title('s(t)');

                % pulse shape plot
                figure; 
                stem(obj.pulseShape);
                title('pulse shape');

                filename = sprintf('bpsk_mod_l%i_upRC%i_p%i_pp%i.dat', length(a_k), N, n_preamble, pos_preamble);
                disp(sprintf('i:filename = %s', filename));
                write_float_binary(st, filename);

            end

        end
        
        function [bitstream] = dem(obj, r_t, T, exp_symb, debug)
            %dem Perform the BPSK demodulation using the given parameters: [bitstream] = dem(obj, r_t, T, exp_symb, debug)
            %
            %   in:
            %
            %   r_t = received sampled signal
            %   T = Symbol period in seconds
            %   exp_symb = # expected symbols
            %   debug = 1 or 0, if 1 plots the first 150 symbols and save a file of st
            %
            %   out:
            %
            %   bitstream = received bitstream after demodulation
            %

            % settings
            delta = T/obj.N; % p(t) sampling period
            p_t = conj(fliplr(obj.pulseShape));

            % correlated with the known received pulse
            % result of the matched filter, length = exp_symb
            ni = zeros(1, exp_symb);

            % length of the raised-cosine = 2*N*GD+1, where N = Upsampling factor,
            % GD = Number of symbol periods between the start of the Raised-cosine
            % and its peak, because it:
            %     i = 0+j*N:((2*N*4+1)+j*N)-1
            for j = 0:length(ni)-1
                for i = 0:(2*obj.N*obj.GD+1)-1
                    ni(j+1) = ni(j+1) + p_t(i+1)*r_t(i+1+j*obj.N);
                end
                ni(j+1) = delta*ni(j+1);
            end
            
            % decision
            
            sign_wave = sign(ni);
            bitstream = zeros(exp_symb, 1);

            for i=1:length(sign_wave)

                if sign_wave(i) == 1
                    bitstream(i) = 1;
                else     
                    bitstream(i) = 0;
                end

            end

            % debug mode
            if debug == 1

                if length(bitstream) > 150
                    toEnd = 150;
                else
                    toEnd = length(bitstream);
                end

                % pulse shape plot
                figure;
                stem(p_t, 'bs');
                title('pulse shape');

                % r(t) plot
                figure;
                stem(r_t(1:(2*obj.N*obj.GD+1)+obj.N*(toEnd-1)));
                title('r(t)');

                % bitstream
                figure;
                stem(bitstream(1:toEnd));
                title('bitstream');

            end
            
        end
        
        function [sym] = bpskMod(obj, a_k, n_preamble, pos_preamble)
            % Perform the bpsk Modulation 0->-1 and 1->1, bits to symbols, [sym] = bpskMod(obj, a_k, n_preamble, pos_preamble)
            % in:
            %   - a_k = Bits to be transmitted
            %   - n_preamble = number of zeros before the marker sequence
            %   - pos_preamble = number of zeros after the marker sequence,
            %   if the n_preamble is zero this value is ignore
            % out:
            %   - sym = bpsk symbol sequence representation of a_k
            
            if size(a_k, 1) > 1
               a_k = a_k'; 
            end
            
            sym = 2*a_k-1;
            
            preamble = [];
            if n_preamble >= 0
                preamble = [zeros(1, n_preamble) (2*obj.know_seq-1)];
                
                if pos_preamble > 0 
                   preamble = [preamble zeros(1, pos_preamble)]; 
                end
            end
            
            sym = [preamble sym];
                
        end
        
        function [st, es] = pulseShaping(obj, sym)
            % perform the pulse shaping operation,  [st, es] = pulseShapping(obj, a_k)
            % in:
            %   - sym = bpsk symbol sequence representation of a_k
            % out:
            %   - st = signal suitable for transmission
            % based on http://gaussianwaves.blogspot.com/2011/04/square-root-raised-cosine-filter.html
            
            sym_os = upsample(sym, obj.N); %oversample sym by N
            
            fprintf('length = %i\n', length(sym_os))
            
            st = conv(sym_os, obj.pulseShape);
%             es = signalEnergy(obj.pulseShape, 1, Inf, length(obj.pulseShape), 0);
            es = pow(obj.pulseShape);
            
        end
        
        function [rtf] = mFilter(obj, r_t)
            % perform the matched filtering. [rtf] = mFilter(r_t)
            % in:
            %   - r_t = noise input signal
            % out:
            %   - rk = matched filtered
            % based on http://gaussianwaves.blogspot.com/2011/04/square-root-raised-cosine-filter.html
            
            rtf = filter(fliplr(obj.pulseShape), 1, r_t);
        end
        
        function [rk] = downSampling(obj, rtf)
            % perform the downsampling to symbols rate at the receiver, [rk] = downSampling(rtf)
            % in:
            %   - rtf = matched filtering noise input signal
            % out:
            %   - rk = downsampled signal
            % based on http://gaussianwaves.blogspot.com/2011/04/square-root-raised-cosine-filter.html
            midSample = length(-obj.GD:1/obj.N:obj.GD);
            %Remove unwanted portions(first few samples till the peak value)
            rtf_truncated = rtf(midSample-1:end); 
            %Now the first sample contains the peak value of the response. From here
            %the samples are extracted depending on the oversampling factor
            rk = downsample(rtf_truncated, obj.N, 1); %here offset=1 means starting from 1st sample %retain every 8th sample
        end
        
%         function [rk] = matchedFiltering(obj, r_t)
%             % perform the matched filtering and downsampling to symbols rate at the receiver, [rk] = matched_filtering(r_t)
%             % in:
%             %   - r_t = noise input signal
%             % out:
%             %   - rk = matched filtered and downsampled signal
%             % based on http://gaussianwaves.blogspot.com/2011/04/square-root-raised-cosine-filter.html
%             
%             y = filter(fliplr(obj.pulseShape), 1, r_t);
%             
%             midSample = length(-obj.GD:1/obj.N:obj.GD);
%             
%             [m, i] = max(y(1:midSample+floor(obj.N/2)));
%             
%             %Remove unwanted portions(first few samples till the peak value)
% %             y_truncated = y(midSample-1:end); 
%             y_truncated = y(i:end);
%             %Now the first sample contains the peak value of the response. From here
%             %the samples are extracted depending on the oversampling factor
% %             rk = downsample(y_truncated, obj.N, 1); %here offset=1 means starting from 1st sample %retain every 8th sample
%             rk = downsample(y_truncated, obj.N);
%         end
        
        function [rk] = matchedFiltering(obj, r_t)
            % perform the matched filtering and downsampling to symbols rate at the receiver, [rk] = matched_filtering(r_t)
            % in:
            %   - r_t = noise input signal
            % out:
            %   - rk = matched filtered and downsampled signal
            % based on http://gaussianwaves.blogspot.com/2011/04/square-root-raised-cosine-filter.html
            
            y = filter(fliplr(obj.pulseShape), 1, r_t);
            
            midSample = length(-obj.GD:1/obj.N:obj.GD);
            %Remove unwanted portions(first few samples till the peak value)
            y_truncated = y(midSample-1:end); 
            %Now the first sample contains the peak value of the response. From here
            %the samples are extracted depending on the oversampling factor
            rk = downsample(y_truncated, obj.N, 1); %here offset=1 means starting from 1st sample %retain every 8th sample
        end
        
        function [bitstream] = decision(obj, rk)
            % Perform the bpsk decision part at the receiver, [bitstream] = decision(rk)
            % in:
            %   - rk = samples obtained after the pulse shape filtering and
            %   downsampling to symbol frequency process perform in
            %   obj.matched_filtering(r_t)
            % out:
            %   - bitstream = decision vector of 1's and 0's, using the
            %   threshold optimum decision
            
            bitstream = (sign(rk)+1)/2;
            
        end
            
    end
    
end

