'''
Copyright (c) 2011, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Feb, 2011
'''

import socket
import time
import sys
import os
from transmitter.TXApp import TXApp
from util.PropertyReader import readProperties

class SlaveTX(object):
    '''
    SlaveTX is responsible of control the TX USRP node.
    '''

    def __init__(self, host, port, path, cr):
        '''
        Constructor
        
        @param host: refers to the local host address
        @param port: port for the server to listen
        @param path: path to the file where the data is stored
        @param cr: path to the file that contains the carrier recovery signal 
        '''
        
        # server
        self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server.bind((host, port))
        self.server.listen(1)
        self.path = path
        self.lc = 0
        
        # carrier recovery
        self.cr = cr
    
    def setTXProperties(self, lo_off, fc, inter, gain, sync, lc):
        '''
        Set the USRP TX properties
        
        @param lo_off: local oscillator offset (int)
        @param fc: Center frequency (float)
        @param inter: Interpolation factor (int)
        @param gain: Gain of the transmitter in dB (int)
        @param sync: True if the Hardware will use the GPSDO (boolean)
        @param lc: value that will be added to the transmitted signal when using Large Carrier method
        '''
        self.lo_off = lo_off
        self.fc = fc
        self.inter = inter
        self.gain = gain
        self.sync = sync
        
        if float(lc) > 0: 
            self.lc = lc
        else:
            self.lc = 0
    
    def launch(self):
        '''
        calls startup for SlaveTx
        '''
        print("i: launch SlaveTX")
        
        while True:
        
            sc, addr = self.server.accept()
            sc.settimeout(10*60)
            
            print("\n\ni: SlaveTX Connection from " + str(addr) + ", time " + time.strftime("%d-%m-%y/%H:%M:%S"))
            
            tic = time.time()
            
            try:
                self.__startup__(sc, addr)
            
            except Exception, e:
                print("e: " + str(e))
            
            sc.close()
            print("i: SlaveTX Connection closed, duration: " + str(time.time() - tic) + " [seg]\n\n")
            
        print("i: SlaveTX end launch")
        
    
    def send(self, sig, isSine, txtime, at):
        '''
        @param sig: filename of the signal (string)
        @param isSine: whether it should transmit a Sine wave or not (bool)
        @param txtime: time of transmission (int)
        @param at: attenuation factor (float)
        
        Responsible for sending data using the USRP
        '''
        
        filename = "" if sig == None else (self.path + sig)
        
        print("i: filename = " + filename)
        
        if sig != None or isSine:
            app = TXApp(self.fc, self.inter, self.gain, "addr=192.168.10.2", at, filename, txtime, self.sync, isSine, self.lc, self.lo_off)
            app.launch()
    
    def __startup__(self, sc, addr):
        '''
        @param sc: socket object
        @param addr: address info 
        
        Responsible for starting the application; for creating and showing
        the initial GUI.
        '''
        print("i: startup")
        
        msg = sc.recv(1024)
        if msg == "start":
            sc.send("ok")
            print("i: start ok")
            
            signal = None
            msg = sc.recv(1024)
            print("i: msg = " + msg)
            while msg != "finish":
                
                tic = time.time()
                
                if msg.find("setSig") >= 0:
                    print("i: setSig received")
                    signal = msg[msg.find(":")+1:msg.rfind(":")]
                    print("i: " + signal)
                    sc.send("ok")
                
                elif msg.find("sendSine") >= 0:
                    print("i: sendSine received")
                    txtime = int(msg[msg.find(":")+1:msg.rfind(":")])
                    self.send(None, True, txtime, 1)
                    sc.send("ok")
                
                elif msg.find("sendCR") >= 0:
                    print("i: sendCR received")
                    self.send(self.cr, False, 0, 1)
                    sc.send("ok")
                
                elif msg.find("sendSignal") >= 0:
                    print("i: sendSignal received")
                    values = msg.split(":")
                    txtime = int(values[1])
                    at = float(values[2])
                    self.send(signal, False, txtime, at)
                    sc.send("ok")
                    
                else:
                    print("i: ending")
                    break

                print("i: cmd duration: " + str(time.time() - tic) + " [seg]\n")

                msg = sc.recv(1024)
                    
        else:
            print("e: not start")
            sc.send("error")
        
        if msg == "finish":
            print("i: finish cmd received")
        
        sc.close()
        print("i: end startup")
        
    def __exit__(self):
        '''
        This method runs on the event dispatching thread.
        '''
        print("somebody call me!")
#        self.__exit__()
        
        
if __name__ == '__main__':
    '''
    Creates an instance of the specified {@code Application}
    subclass, sets the {@code ApplicationContext} {@code
    application} property, and then calls the new {@code
    Application's} {@code startup} method.  The {@code launch} method is
    typically called from the Application's {@code main}:
    '''
        
    # Reading the properties
    confFile = "confTX.txt"
    if(len(sys.argv) > 1):
        arg = sys.argv[1]
        confFile = arg if len(arg) > 0 else confFile 
    else:
        print("working with default config file path")
    
    properties = readProperties(confFile)
    
    print("Properties:")
    for p in properties:
        print("\t" + p + " : " + properties[p])
    
    # fixing/checking the properties   
    path = properties["txpath"] 
    path = path if (path.endswith("/")) else path+"/" 
    
    if not os.path.exists(path + properties["crsig"]):
        properties["crsig"] = None
        print("i: unable to find the Carrier recovery signal. Nothing will be send instead")
    
    sync = True if properties["sync"]=="True" else False
    
    # creating the app
    app = SlaveTX(properties["txip"],
                  int(properties["txport"]),
                  path,
                  properties["crsig"])
    
    app.setTXProperties(int(properties["lo_off"]), 
                        float(properties["fc"]), 
                        int(properties["inter"]),
                        int(properties["gain"]), 
                        sync,
                        float(properties["lc"]))
    
    app.launch()
    exit()