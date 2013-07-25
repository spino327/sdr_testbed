'''
Copyright (c) 2010, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Dec, 2010
'''

from gnuradio import gr
from gnuradio import eng_notation

n2s = eng_notation.num_to_str

class RXMatchedFilter(gr.hier_block2):
    '''
    This class performs the matched filer (signal demodulator) for reception
    
    
        --->(rrc)--->

    '''

    def __init__(self, gain, samp_freq, sym_rat, roll_off, pern):
        '''
        in:
            - gain = gain of the filter
            - samp_freq = sampling rate of each pulse period, (samples per symbol)
            - sym_rat = per n times samp freq, how much symbols?
            - roll_off = squared root raised cosine's roll off factor
            - pern = periods number
        '''
        gr.hier_block2.__init__(self, "Pulse shaping",
                                gr.io_signature(1, 1, gr.sizeof_gr_complex),
                                gr.io_signature(1, 1, gr.sizeof_gr_complex))
        
        # instance variables
        sps = int(samp_freq/sym_rat)
        print "i: IF sample rate:", n2s(samp_freq)
        print "i: Symbol rate:", n2s(sym_rat)
        print "i: Samples/symbol:", sps
        print "i: RRC bandwidth:", roll_off
        
        # matched filter
        taps = gr.firdes.root_raised_cosine(gain, samp_freq, sym_rat, roll_off, pern*samp_freq)
        self.rrc = gr.interp_fir_filter_ccf(1, # Interpolation rate
                                            taps)        
        
        # connections
        self.connect((self, 0), (self.rrc, 0), (self, 0))
        
        