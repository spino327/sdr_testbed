'''
Copyright (c) 2011, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Apr, 2011
'''

from gnuradio import gr
from gnuradio.gr import firdes
from gnuradio import blks2
from util import USRP2Conf

class TxUSRP(gr.hier_block2):
    '''
    
    This class handle the samples rate fixing operation and also the frequency error fixing operation
    
    Resampler a lower signal rate to the requirement rate at the usrp
    
        --->(pfb_resampler)--->(xlatting_filter)--->(usrp_sink)
    
    '''

    def __init__(self, *params):
        gr.hier_block2.__init__(self, "TxUSPR",
                                gr.io_signature(1, 1, gr.sizeof_gr_complex),
                                gr.io_signature(0, 0, 0))
        
        if len(params) == 7:
            self.__uhd(params[0], params[1], params[2], params[3], params[4], params[5], params[6])
#        elif len(params) == 6:
#            self.__raw(params[0], params[1], params[2], params[3], params[4], params[5])
        else:
            raise Exception()


    def __uhd(self, fc, lo_off, inter, gain, addr, inSampRate, sync):
        '''
            in: 
                - fc = center frequency
                - lo_off = LO off
                - inter = interporlation factor
                - gain = gain in the tx, only with 2450
                - addr = ip address, format = "addr=ip, mimo_mode="
                - inSampRate = incoming sample frequency, basically here we determine the re-sampler interpolation factor 
                - sync = True is we're gonna use an external ref clock
        
        '''
        
        # instance variables
        self.type = "UHD"
        (self.tx, basebandFreq, dxcFreq) = USRP2Conf.getUhdUSRPSink(fc, lo_off, inter, gain, addr, sync)
        sampRate = float(self.tx.get_clock_rate())/inter
        
        self.postProcessing(inSampRate, dxcFreq, sampRate)


#    def __raw(self, fc, inter, gain, eth, inSampRate, sync):
#        '''
#            in: 
#                - fc = center frequency
#                - inter = interporlation factor
#                - gain = gain in the tx, only with 2450
#                - eth = ethernet interface name(String)
#                - inSampRate = incoming sample frequency, basically here we determine the re-sampler interpolation factor 
#                - sync = True is we're gonna use an external ref clock
#        
#        '''
#        
#        # instance variables
#        self.type = "RAW"
#        (self.tx, basebandFreq, dxcFreq) = USRP2Conf.getUSRP2Sink(fc, inter, gain, eth, sync)
#        sampRate = float(self.tx.dac_rate())/inter
        
#        self.postProcessing(inSampRate, dxcFreq, sampRate)

    def postProcessing(self, inSampRate, dxcFreq, sampRate):
        
        # xlating
        if dxcFreq != 0:
            xlateFilterTaps = firdes.low_pass(1, sampRate, sampRate / 2, sampRate / 10, firdes.WIN_HAMMING, 6.76)
            self.xlatingFilter = gr.freq_xlating_fir_filter_ccc(1, (xlateFilterTaps), 
                dxcFreq, 
                sampRate)
            print "i: xlating filter fixed to " + str(dxcFreq)
        else:
            self.xlatingFilter = gr.multiply_const_vcc((1, ))
            print "i: xlating filter not needed"
            
        # pfb resampler
        self.resamplerFactor = sampRate / inSampRate
        nphases = 32
        frac_bw = 0.45
        rs_taps = firdes.low_pass(nphases, nphases, frac_bw, 0.5 - frac_bw)
        self.resampler = blks2.pfb_arb_resampler_ccf(self.resamplerFactor, 
            (rs_taps), 
            nphases)
        print "i: re-sampler relation new_freq/old_freq = " + str(self.resamplerFactor)
        #EO instance variables
        
        self.isRTEnable = gr.enable_realtime_scheduling()
        if self.isRTEnable == gr.RT_OK:
            print "i: realtime enable: True"
        else:
            print "i: realtime enable: False"
        
        # Connections
        self.connect((self, 0), (self.resampler, 0), (self.xlatingFilter, 0), (self.tx, 0))

    def dac_rate(self):
        '''
        return the DAC rate in Hz
        '''
        if self.type == "UHD":
            return self.tx.get_clock_rate()
        else:
            return self.tx.dac_rate()
        
        