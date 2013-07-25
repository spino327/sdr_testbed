'''
Copyright (c) 2010, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Nov, 2010
'''

from util.RXMatchedFilter import RXMatchedFilter

from gnuradio import gr
from gnuradio import blks2
from gnuradio.gr import firdes


class RXData(gr.hier_block2):
    '''
    @fixme must update the processing blocks diagram
    
    RXData represents the data processing general block
    
    --->(bandpass)--->(resampler)--->(FileSink)
    
    '''

    def __init__(self, app, samp_rate, inter, dec, filename):
        '''
        in:
            - app = object of type RXApp
            - samp_rate = sample rate in Hertz
            - inter = interpolation factor
            - dec = decimation factor
            - filename
            
        output
            - port 0 = MF
        '''
        gr.hier_block2.__init__(self, "RXData",
                                gr.io_signature(1, 1, gr.sizeof_gr_complex),
                                gr.io_signature(0, 0, 0))
        
        # instance variables
        self.filename = filename
        self.isRecording = False
        self.app = app
        
        #  "bandpass" to new sample rate
        # cleaning the out-of-band interferences, low pass filter around the shape 
        lowPassTaps = firdes.low_pass(1, samp_rate, samp_rate/2, samp_rate/20, firdes.WIN_HAMMING, 6.76)
        self.lowPass = gr.multiply_const_vcc([1, ]) #gr.fir_filter_ccc(1, lowPassTaps)
        
        #Matched Filter
#        self.matched = gr.multiply_const_cc(1)   #RXMatchedFilter(1, 20, 1, 0.8, 10)
        
        # resampler after Matched-filter
        self.resampler = blks2.rational_resampler_ccc(
                interpolation=inter,
                decimation=dec,
                taps=None,
                fractional_bw=None,
        )
        
        # saving the IQ samples
        self.fileSink = gr.file_sink(gr.sizeof_gr_complex*1, self.filename) if(self.filename <> None and len(self.filename) > 0) else None
        
        #EO instance variables
        
        self.__makeConnections()

        
    def __makeConnections(self):
        '''
        uses the method connect(src, des)
        '''

        # smart receiver
        self.connect((self, 0), (self.lowPass, 0))
#        self.connect((self.lowPass, 0), (self.matched, 0))
#        self.connect((self.matched, 0), (self.resampler, 0))
        self.connect((self.lowPass, 0), (self.resampler, 0))
        # mf samples
        self.connect((self.resampler, 0), (self.fileSink, 0))
        
        