'''
Copyright (c) 2011, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Jul, 2011
'''

from receiver.RxUsrp import RxUSRP
from receiver.RXData import RXData
from receiver.RXControlApp import RXControlApp
from gnuradio import gr
import time

class RXApp(object):
    '''
    Application life cycle
    '''
    
    def __init__(self, fc, dec, gain, eth, sync, filename, *params):
        '''
        Create the instances of the gui and data:
        in:
            - fc = center frequency at bandpass
            - dec = usrp decimation factor
            - gain = RF gain
            - eth = name of the ethernet
            - sync = True if u wanna sync with external clock
        '''
        # desired signal sample rate
        outSampRate = float(100e6)/dec
        
        if len(params) > 0:
            self.rx = RxUSRP(fc, params[0], dec, gain, eth, outSampRate, sync)
            
        self.data = RXData(self, outSampRate, 1, 1, filename)
        self.tb = gr.top_block()
        
#        if recording_time <= 0:
#            recording_time = 1
#        
#        self.control = RXControlApp(self, recording_time) if(recording_time > 0) else None
        
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
        
        # from the antenna to the file
        self.tb.connect((self.rx, 0), (self.data, 0))
        
#        if self.control <> None:
#            self.control.start()
        
        # running
        self.tb.start()
        
        print "i: end startup"
        
    def stopApp(self):
        """
        Stop app
        """
        try:
            self.tb.stop()
            
        except Exception, e:
            print("e: ", e)
            return False
            
        return True
