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
from transmitter.TXData2 import TXData2
from transmitter.TxUsrp import TxUSRP
from transmitter.TXControlApp import TXControlApp
from gnuradio import gr

class TXApp(object):
    '''
    Application life cycle
    '''

    def __init__(self, fc, inter, gain, eth, at, filename, txtime, sync, sine, lc, *params):
        '''
        Create the instances of the gui and data:
        in:
            - fc = center frequency at bandpass
            - inter = usrp interpolation factor
            - gain = RF gain
            - eth = name of the ethernet
            - at = attenuation factor
            - filename = filename
            - txtime = in seconds, txtime=0 if you want to send it once, txtime>0 if you want to send it for x sec  
            - sync = True if u wanna sync with external clock
            - sine
            - params = self.lo_off
        '''
        
        # desired signal sample rate
        inSampRate = float(100e6)/inter #44.1e3
        
        if len(params) > 0:
            self.tx = TxUSRP(fc, params[0], inter, gain, eth, inSampRate, sync)
        
        repeat = True if(txtime > 0) else False

        self.data = TXData2(self, self.tx.dac_rate()/inter, at, filename, repeat, sine, lc)
        self.tb = gr.top_block()
        
        self.control = TXControlApp(self, txtime) if(txtime > 0) else None
        
    def launch(self):
        '''
        calls startup
        '''
        print "i: launch"
        self.__startup__()
        print "i: end launch"
    
    def __startup__(self):
        '''
        Responsible for starting the application; for creating and showing
        the initial GUI.
        '''
        print "i: startup"
        
        # from the file to the antenna 
        self.tb.connect((self.data, 0), (self.tx, 0))
        
        if self.control <> None:
            self.control.start()
        
        # running
        self.tb.run()
            
        print "i: end startup"
        
