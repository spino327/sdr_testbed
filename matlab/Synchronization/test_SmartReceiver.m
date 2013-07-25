classdef test_SmartReceiver < mlunit.test_case
    %TESTSMARTRECEIVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        verbose = true
    end
    
    methods
        
        % constructor
        function self = test_SmartReceiver(varargin)
            self = self@mlunit.test_case(varargin{:})
            
        end
        
        function self = testOldModFindMarkerSequenceZeroNoise(self)

            rcosine = importdata('data/rcosine/srrc_up20_gd5_ro20.mat');
            a_k = importdata('data/a_k/a_k72.mat');
            nmsg = 5; % number of msg repetitions
            expect = nmsg;
            ks = [1 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 1 1 1];%[1 1 1 1 1 1 1 1 1 1 0 0];
            preamble = 100; % number of zero value signal samples
            N = 20; % upsampling factor
            
            bpskObj = ModUtils(rcosine, N, ks);
            s_t = bpskObj.mod(a_k, preamble, 0, 0)';
            header = bpskObj.mod(ks, 0, 0, 0)';
            
            data = [];
            for i=1:nmsg
                data = [data s_t];
            end
            
            % finding marker sequences
            smRec = SmartReceiver(header);
            headstart = smRec.findMarkerSequence(data, 1);
            result = length(headstart);
            mlunit.assert_equals(expect, result);
            disp('i:finding marker sequences - OK')
            
            % comparing the position of the marker signals
            expect = sort(preamble*N + [0:nmsg-1]*length(s_t));
            result = sort(headstart);
            mlunit.assert_equals(expect, result);
            disp('i:comparing the position of the marker signals - OK')
        end
        
        function self = testOldModFindMarkerSequenceNoisy(self)

            rcosine = importdata('data/rcosine/rrc_up20_gd5_GNU2.mat');
            a_k = importdata('data/a_k/rand_bits_1E3.mat');
            nmsg = 5; % number of msg repetitions
            expect = nmsg;
            ks = [1 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 1 1 1];%[1 1 1 1 1 1 1 1 1 1 0 0];
            preamble = 100; % number of zero value signal samples
            N = 20; % upsampling factor
            
            bpskObj = BPSK(rcosine, N, ks);
            s_t = bpskObj.mod(a_k, preamble, 0, 0)';
            header = bpskObj.mod(ks, 0, 0, 0)';
     
            data = [];
            for i=1:nmsg
                data = [data s_t];
            end
            
            vx = 0 + (0.1-0).*rand();
            data = data + sqrt(vx)*randn(size(data));
            disp(sprintf('i: vx = %f', vx));
            
            % finding marker sequences
            smRec = SmartReceiver(header);
            [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence(data, 1);
            
            if self.verbose
                figure('Name', 'testOldModFindMarkerSequenceNoisy');
                subplot(2,1, 1); hold on; plot(data); stem(headstart, data(headstart), 'rs'); hold off; title(sprintf('s(t), vx=%f', vx));
                subplot(2,1, 2); plot(ycorr); title('xcorr')
            end
            
            result = length(headstart);
            mlunit.assert_equals(expect, result);
            disp('i: finding marker sequences - OK')
            
            % comparing the position of the marker signals
            expect = sort(preamble*N + [0:nmsg-1]*length(s_t));
            result = sort(headstart);
            % using an error tolerance, less than 5 samples
            error = norm(expect-result);
            mlunit.assert(error<nmsg*5);
            disp('i:comparing the position of the marker signals - OK')
        end
        
        function self = testNewModFindMarkerSequenceZeroNoise(self)

            rcosine = importdata('data/rcosine/rrc_up20_gd5_GNU2.mat');
            a_k = importdata('data/a_k/rand_bits_1E3.mat');
            nmsg = 5; % number of msg repetitions
            expect = nmsg;
            ks = [1 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 1 1 1];%[1 1 1 1 1 1 1 1 1 1 0 0];
            preamble = 100; % number of zero value signal samples
            N = 20; % upsampling factor
            
            bpskObj = BPSK(rcosine, N, ks);
            % s(t)
            [sym] = bpskObj.bpskMod(a_k, preamble, 0);
            s_t = bpskObj.pulseShaping(sym);
            % pn
            [sym] = bpskObj.bpskMod(ks, 0, 0);
            header = bpskObj.pulseShaping(sym);
            
            data = [];
            for i=1:nmsg
                data = [data s_t];
            end
            
            % finding marker sequences
            smRec = SmartReceiver(header);
            headstart = smRec.findMarkerSequence(data, 1);
            result = length(headstart);
            mlunit.assert_equals(expect, result);
            disp('i:finding marker sequences - OK')
            
            % comparing the position of the marker signals
            expect = sort(preamble*N + [0:nmsg-1]*length(s_t));
            result = sort(headstart);
            mlunit.assert_equals(expect, result);
            disp('i:comparing the position of the marker signals - OK')
        end
        
        function self = testNewModFindMarkerSequenceNoisy(self)

            rcosine = importdata('data/rcosine/rrc_up20_gd5_GNU2.mat');
            a_k = importdata('data/a_k/rand_bits_1E3.mat');
            nmsg = 5; % number of msg repetitions
            expect = nmsg;
            ks = [1 1 1 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 1 1 1];
            preamble = 100; % number of zero value signal samples
            N = 20; % upsampling factor
            
            bpskObj = BPSK(rcosine, N, ks);
            % s(t)
            [sym] = bpskObj.bpskMod(a_k, preamble, 0);
            s_t = bpskObj.pulseShaping(sym);
            % pn
            [sym] = bpskObj.bpskMod(ks, 0, 0);
            header = bpskObj.pulseShaping(sym);
     
            data = [];
            for i=1:nmsg
                data = [data s_t];
            end
            vx = 0 + (0.1-0).*rand();
            data = data + sqrt(vx)*randn(size(data));
            disp(sprintf('i: vx = %f', vx));
            % finding marker sequences
            smRec = SmartReceiver(header);
            [headstart, ycorr, peaks, loc] = smRec.findMarkerSequence(data, 1);
            
            if self.verbose
                figure('Name', 'testNewModFindMarkerSequenceNoisy'); 
                subplot(2,1, 1); hold on; plot(data); stem(headstart, data(headstart), 'rs'); title('s(t)'); hold off;
                subplot(2,1, 2); plot(ycorr); title('xcorr')
            end
            
            result = length(headstart);
            mlunit.assert_equals(expect, result);
            disp('i:finding marker sequences - OK')
            
            % comparing the position of the marker signals
            expect = sort(preamble*N + [0:nmsg-1]*length(s_t));
            result = sort(headstart);
            % using an error tolerance, ber
            [num, ratio] = biterr(expect, result);
            disp(sprintf('i: nerr = %i, ber = %f', num, ratio));
            mlunit.assert(ratio < 0.2);
            disp('i:comparing the position of the marker signals - OK')
        end
        
        function self = testFreqSync(self)
            
            st = read_float_binary('data/carrier_recovery/bpsksc.dat');
            
            phoff = 2*pi*rand();
            
            rt = st.*exp(1i*phoff);

            objSR = SmartReceiver([]);
            rtFix = objSR.frequencySync(rt);
            
            if self.verbose
                figure('Name', 'testFreqSync-Noiseless');
                subplot(3, 1, 1); plot(st, 'r'); title('s(t)');
                subplot(3, 1, 2); hold on; plot(real(rt), 'r'); plot(imag(rt), 'b'); title('r(t)'); hold off;
                subplot(3, 1, 3); hold on; plot(real(rtFix), 'r'); plot(imag(rtFix), 'b'); title('rtFix'); hold off;
            end
            
            mlunit.assert(1-pow(real(rtFix))/pow(real(st))<0.1);
            disp('i:freq sync zero noise - OK')
            
            vx = 0 + (0.1-0).*rand();
            rt = rt + sqrt(vx)*randn(1, length(rt))';
            
            rtFix = objSR.frequencySync(rt);
            
            if self.verbose
                figure('Name', 'testFreqSync-Noise');
                subplot(3, 1, 1); plot(st, 'r'); title('s(t)');
                subplot(3, 1, 2); hold on; plot(real(rt), 'r'); plot(imag(rt), 'b'); title('r(t)'); hold off;
                subplot(3, 1, 3); hold on; plot(real(rtFix), 'r'); plot(imag(rtFix), 'b'); title('rtFix'); hold off;
            end
            
            mlunit.assert(1-pow(real(rtFix))/pow(real(st))<0.1);
            disp('i:freq sync noise - OK')
            
        end
        
    end
    
end

