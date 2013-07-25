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

from gnuradio import gr
from gnuradio import blks2
from gnuradio.gr import firdes

class TXData(gr.hier_block2):
    '''
    TXData represents the data processing general block
    
        (filesrc)--->(multiply)--->(f2c)--->
   
    '''

    def __init__(self, app, samp_rate, at, filename, repeat, sine):
        '''
        in:
            - app = object of type RXApp
            - samp_rate = sample rate in Hertz
            - at = attenuation factor
            - filename = filename
            - repeat = if True them reads in a loop the file
            - sine
        '''
        gr.hier_block2.__init__(self, "RXData",
                                gr.io_signature(0, 0, 0),
                                gr.io_signature(1, 1, gr.sizeof_gr_complex))
        
        
        # instance variables
        self.app = app
        
        if sine:
            self.fileSrc = gr.sig_source_f(samp_rate, gr.GR_SIN_WAVE, 1000, 1.0)
        else:
            self.fileSrc = gr.file_source(gr.sizeof_float*1, filename, repeat)
#        self.fileSrc = gr.sig_source_f(samp_rate, gr.GR_SIN_WAVE, 1000, 1.0)
#        self.fileSrc = gr.sig_source_f(samp_rate, gr.GR_CONST_WAVE, 1000, 1.0)
        self.mulitplyCte = gr.multiply_const_vff((at, ))
        self.f2c = gr.float_to_complex(1)
        
        #EO instance variables
        
        self.__makeConnections()
        
        
    def __makeConnections(self):
        '''
        uses the method connect(src, des)
        '''
        self.connect((self.fileSrc), (self.mulitplyCte, 0))
        self.connect((self.mulitplyCte, 0), (self.f2c, 0))
        self.connect((self.f2c, 0), (self, 0))
        