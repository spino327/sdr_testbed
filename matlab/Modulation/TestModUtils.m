%Copyright (c) 2011, Universidad Industrial de Santander, Colombia
%University of Delaware
%All rights reserved.
%@author: Sergio Pino
%@author: Henry Arguello
%Website: http://www.eecis.udel.edu/
%emails  : sergiop@udel.edu - henarfu@udel.edu
%Date   : Jun, 2011

classdef TestModUtils < mlunit.test_case
    %TESTMODUTILS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        verbose = false
        verboseSDR = false
        verboseSDR2 = false
        isRecording = true
    end
    
    methods
        
        function self = TestModUtils(varargin)
            self = self@mlunit.test_case(varargin{:});
            
        end
        
        function self = testPulseShape(self)
            
            tic
            a_k = importdata('data/a_k/a_k1E4.mat');%[1 0 0 0 0 1 1 1 1 1 0 1 0];%importdata('data/a_k/a_k72.mat');%
            ks = [];
            
            N = 20; % upsampling factor
            gd = 5;
            ro = 20;
            
            ps = importdata(sprintf('data/rcosine/srrc_up%i_gd%i_ro%i.mat', N, gd, ro)); % pulse shape
            
            modUtilsObj = ModUtils(ps, N, []);
            
            expected = modUtilsObj.mod(a_k, 0, 0, 0);
            
            [sym] = modUtilsObj.bpskMod(a_k, 0, 0);
            [result, es] = modUtilsObj.pulseShaping(sym);
            fprintf('es = %f\n', es);
            fprintf('result length = %f\n', length(result));
            
            toc
            if self.verbose
                figure('Name', 'testPulseShape'); 
                hold on; stem(expected, 'r'); plot(result, 'b'); 
                hold off; legend('expected', 'result'); title('Pulse Shaping');
            end
            
            mlunit.assert_equals(expected, result(1:length(expected))');
            disp('i: pass pulse shaping')
                        
            if self.isRecording
               save(sprintf('data/Testdata/testPulseShape%i.mat', length(a_k)), 'result', '-v7');
               disp('i: pass pulse shaping')
            end
            
        end
        
        function self = testNewModDemTransmissionZeroNoise(self)
            
            tic
            % creating a sample file of the modulated signal
            N = 20; % upsampling factor
            gd = 5;
            ro = 20;
            folder = 'data/tx/12-02-11/';

            ps = importdata(sprintf('data/rcosine/srrc_up%i_gd%i_ro%i.mat', N, gd, ro)); % pulse shape
            ps = ps/sqrt(sum((ps).^2));


            a_k = importdata('data/a_k/a_k1E4.mat');
            ks = importdata('data/ks/ks101.mat');

            n_preamble = 500;
            pos_preamble = 0;

            % mod
            bpskObj = ModUtils(ps, N, ks);
            [sym] = bpskObj.bpskMod(a_k, n_preamble, pos_preamble);
            [st, es] = bpskObj.pulseShaping(sym);
            
            % noise free, matched filtering and decision
            rk = bpskObj.matchedFiltering(st);
            result = bpskObj.decision(rk(n_preamble+1:end));
            
            toc
            if self.verbose
                figure('Name', 'testNewModDemTransmissionZeroNoise'); 
                subplot(3, 1, 1); plot(st); title('s(t)');
                subplot(3, 1, 2); stem(rk); title('r_k');
                subplot(3, 1, 3); hold on; stem(a_k, 'sr'); stem(result); 
                hold off; legend('a_k', 'a_krec'); title('data')
            end
            
            
            mlunit.assert_equals(ks, result(1:length(ks)));
            mlunit.assert_equals(a_k, result(length(ks)+1:end));
            disp('i: received sequence is equal in New ModDem')
            
            if self.isRecording
                % saving
                filename = strcat(folder, sprintf('bpsk_mod_l%i_upRC%i_gd%i_ro%i_p%i_pp%i.dat', ...
                                                    length(a_k), N, gd, ro, n_preamble, pos_preamble));
                fprintf('i:filename = %s\n', filename);
                write_float_binary(st, filename); 
            end
            
        end
        
        function self = testMatchedFiltering2Steps(self)
            
            tic
            % creating a sample file of the modulated signal
            N = 20; % upsampling factor
            gd = 5;
            ro = 20;
            folder = 'data/tx/12-02-11/';

            ps = importdata(sprintf('data/rcosine/srrc_up%i_gd%i_ro%i.mat', N, gd, ro)); % pulse shape
            ps = ps/sqrt(sum((ps).^2));


            a_k = importdata('data/a_k/a_k1E4.mat');
            ks = importdata('data/ks/ks101.mat');

            n_preamble = 500;
            pos_preamble = 0;

            % mod
            bpskObj = ModUtils(ps, N, ks);
            [sym] = bpskObj.bpskMod(a_k, n_preamble, pos_preamble);
            [st, es] = bpskObj.pulseShaping(sym);
            
            % noise free, matched filtering and decision
            rk = bpskObj.matchedFiltering(st);
            expected = bpskObj.decision(rk(n_preamble+1:end));
            
            mlunit.assert_equals(ks, expected(1:length(ks)));
            mlunit.assert_equals(a_k, expected(length(ks)+1:end));
            disp('i: received sequence is equal in New ModDem - Normal')
            
            % noise free, matched filtering 2 steps
            rtf = bpskObj.mFilter(st);
            rk = bpskObj.downSampling(rtf);
            result = bpskObj.decision(rk(n_preamble+1:end));
            
            mlunit.assert_equals(ks, result(1:length(ks)));
            mlunit.assert_equals(a_k, result(length(ks)+1:end));
            disp('i: received sequence is equal in New ModDem - 2 Steps')
            
            toc
%             if self.verbose
%                 figure('Name', 'testNewModDemTransmissionZeroNoise'); 
%                 subplot(3, 1, 1); plot(st); title('s(t)');
%                 subplot(3, 1, 2); stem(rk); title('r_k');
%                 subplot(3, 1, 3); hold on; stem(a_k, 'sr'); stem(result); 
%                 hold off; legend('a_k', 'a_krec'); title('data')
%             end
                  
%             if self.isRecording
%                 % saving
%                 filename = strcat(folder, sprintf('bpsk_mod_l%i_upRC%i_gd%i_ro%i_p%i_pp%i.dat', ...
%                                                     length(a_k), N, gd, ro, n_preamble, pos_preamble));
%                 fprintf('i:filename = %s\n', filename);
%                 write_float_binary(st, filename); 
%             end
            
        end
        
    end
    
end

