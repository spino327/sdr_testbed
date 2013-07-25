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
import unittest
from transmitter.TxUsrp import TxUSRP
from gnuradio import gr

class Test_TxUsrp(unittest.TestCase):


    def testCreation(self):
        
        type = "UHD"
        inter = 512
        inSampRate = 50e3
        expResam = (100e6/inter)/inSampRate
        sync = True
        
        if type == "UHD":
            tx = TxUSRP(2400e6, 5000, inter, 0, "addr=192.168.10.2", inSampRate, sync)
        else:
            tx = TxUSRP(2400e6, inter, 0, "eth1", inSampRate, sync)
        
        resResam = tx.resamplerFactor
        self.assertEquals(expResam, resResam)
        
        # Realtime
        expected = gr.RT_OK
        result = tx.isRTEnable
        self.assertEqual(expected, result)
        
        # dac
        expected = 100e6
        result = tx.dac_rate()
        self.assertEqual(expected, result)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()