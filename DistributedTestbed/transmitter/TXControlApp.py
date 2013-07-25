'''
Copyright (c) 2010, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Nov, 2011
'''

import threading
import time

class TXControlApp(threading.Thread):
    '''
    This class control the running of an gnuradio class for TX Slave
    '''
    def __init__(self, app, time):
        '''
        in:
            - app = instance of gr.top_block or a subclass
            - time = desired time in seconds
            - gain = tuple of gains
        '''
        threading.Thread.__init__(self)
        self.app = app
        self.time = time

    def run(self):
        
        print("i: TX running control, waiting " + str(self.time))
        time.sleep(self.time)
        self.app.tb.stop()
        print("i: TX Done, tb shutdown")
