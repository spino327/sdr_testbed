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

import threading
import time

class ControlApp(threading.Thread):
    '''
    This class control the running of an gnuradio class
    '''
    def __init__(self, app, time, gain):
        '''
        in:
            - app = instance of gr.top_block or a subclass
            - time = desired time in seconds
            - gain = tuple of gains
        '''
        threading.Thread.__init__(self)
        self.app = app
        self.time = time
        self.gain = gain

    def run(self):
        #1
#        print("i: running control, waiting " + str(self.time))
#        time.sleep(self.time)
#        self.tb.stop()
#        print("i: Done, tb shutdown")
        
        #2
        self.app.tb.stop()
        self.app.tb.wait()
        self.app.data.stopRecording()
        
        for g in self.gain:
            print "i: recording at gain = " + str(g) 
            self.app.setGain(g)
            time.sleep(1)
            self.app.data.startRecording(g)
            self.app.tb.start()
            print("i: running control, waiting " + str(self.time))
            time.sleep(self.time)
            self.app.tb.stop()
            self.app.tb.wait()
            self.app.data.stopRecording()

        #3
#        print("i: running control, waiting " + str(self.time))
#        time.sleep(self.time)
#        self.app.tb.stop()
#        print("i: Done, tb shutdown")
