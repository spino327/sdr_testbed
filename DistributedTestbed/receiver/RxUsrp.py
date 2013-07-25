'''
Copyright (c) 2010, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Apr, 2010
'''

from gnuradio import gr
from gnuradio.gr import firdes
from gnuradio import blks2
from gnuradio import uhd
from util import USRP2Conf

class RxUSRP(gr.hier_block2):
    '''
    
    This class handle the samples rate fixing operation and also the frequency error fixing operation
    
    Resampler a higher signal rate to the requirement rate at the receiver
    
        (usrp_source)--->(xlatting_filter)--->(pfb_resampler)--->
    
    '''
    
    def __init__(self, *params):
        
        gr.hier_block2.__init__(self, "RxUSPR",
                        gr.io_signature(0, 0, 0),
                        gr.io_signature(1, 1, gr.sizeof_gr_complex))
        
        if len(params) == 7:
            self.__uhd(params[0], params[1], params[2], params[3], params[4], params[5], params[6])
#        elif len(params) == 6:
#            self.__raw(params[0], params[1], params[2], params[3], params[4], params[5])
        else:
            raise Exception()

    def __uhd(self, fc, lo_off, dec, gain, addr, outSampRate, sync):
        '''
            For UHD fw and fpga
            in: 
                - fc = center frequency
                - lo_off = LO off
                - dec = decimation factor
                - gain = gain in the rx, only with 2450
                - addr = ip address, format = "addr=ip, mimo_mode="
                - outSampRate = outgoing sample frequency, basically here we determine the re-sampler decimation/interpolation factor
                - sync = True is we're gonna use an external ref clock
        
        '''
        
        # instance variables
        self.type = "UHD"
        # usrp
        (self.rx, basebandFreq, dxcFreq) = USRP2Conf.getUhdUSRP2Source(fc, lo_off, dec, gain, addr, sync)
        sampRate = float(self.rx.get_clock_rate())/dec
        
        self.postProcessing(outSampRate, dxcFreq, sampRate)
        
        
#   def __raw(self, fc, dec, gain, eth, outSampRate, sync):
#        '''
#            For war ethernet fw and fpga
#            in: 
#                - fc = center frequency
#                - dec = decimation factor
#                - gain = gain in the rx, only with 2450
#                - eth = ethernet interface name(String)
#                - outSampRate = outgoing sample frequency, basically here we determine the re-sampler decimation/interpolation factor
#                - sync = True is we're gonna use an external ref clock
#        
#        '''
#        
#        # instance variables
#        self.type = "RAW"
#        # usrp
#        (self.rx, basebandFreq, dxcFreq) = USRP2Conf.getUSRP2Source(fc, dec, gain, eth, sync)
#        sampRate = float(self.rx.adc_rate())/dec
#
#        self.postProcessing(outSampRate, dxcFreq, sampRate)

    def postProcessing(self, outSampRate, dxcFreq, sampRate):
        
        stXlate = False
        
        # xlating
        if dxcFreq != 0:
            xlateFilterTaps = firdes.low_pass(1, sampRate, sampRate / 2, sampRate / 10, firdes.WIN_HAMMING, 6.76)
            self.xlatingFilter = gr.freq_xlating_fir_filter_ccc(1, (xlateFilterTaps), 
                dxcFreq, 
                sampRate)
            print "i: xlating filter fixed to " + str(dxcFreq)
            
            stXlate = True
            
        else:
#            self.xlatingFilter = gr.multiply_const_vcc((1, ))
            stXlate = False
            print "i: xlating filter not needed"
        
        # pfb resampler
        self.resamplerFactor = outSampRate / sampRate
        
        if self.resamplerFactor != 1:
        
            nphases = 32
            frac_bw = 0.45
            rs_taps = firdes.low_pass(nphases, nphases, frac_bw, 0.5 - frac_bw)
            self.resampler = blks2.pfb_arb_resampler_ccf(self.resamplerFactor, 
                (rs_taps), 
                nphases)
            print "i: re-sampler relation new_freq/old_freq = " + str(self.resamplerFactor)
            #EO instance variables
            
            # Connections
            if stXlate:
                self.connect((self.rx, 0), (self.resampler, 0), (self.xlatingFilter, 0), (self, 0))
            else:
                self.connect((self.rx, 0), (self.resampler, 0), (self, 0))
            
        else:
#            print "i: stXlate = " + str(self.resamplerFactor)
            # Connections
            if stXlate:
                self.connect((self.rx, 0), (self.xlatingFilter, 0), (self, 0))
            else:
                self.connect((self.rx, 0), (self, 0))
        
        #EO instance variables
        
        self.isRTEnable = gr.enable_realtime_scheduling()
        if self.isRTEnable == gr.RT_OK:
            print "i: realtime enable: True"
        else:
            print "i: realtime enable: False"
      
      
    def adc_rate(self):
        '''
        return the ADC rate in Hz
        '''
        if self.type == "UHD":
            return self.rx.get_clock_rate()
        else:
            return self.rx.adc_rate()
        
    
    def set_gain(self, gain):
        '''
        set the gain of the usrp
        '''
        self.rx.set_gain(gain)
        return self.rx.get_gain()
    
    def get_gain(self):
        '''
        get the gain of the usrp
        '''
        return self.rx.get_gain()
        
    def set_center_freq(self, fc, lo_off):
        '''
        set the center frequency
        '''
        tune_request = uhd.tune_request_t(fc, lo_off)
        tune_res = self.rx.set_center_freq(tune_request)
        rf_freq = tune_res.actual_rf_freq
        dxc_freq = 0
        print("i: actual_rf_freq = " + str(rf_freq) + 
              "\ni: actual_dsp_freq = " + str(tune_res.actual_dsp_freq) + 
              "\ni: target_rf_freq = " + str(tune_res.target_rf_freq) +
              "\ni: target_dsp_freq = " + str(tune_res.target_dsp_freq))
        
#        return (rf_sink, rf_freq, dxc_freq)

        