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
from receiver.RXApp import RXApp
from util.PropertyReader import readProperties
from util import Utils

class SlaveRX(object):
    '''
    SlaveRX is responsible of control the RX USRP node.
    '''

    def __init__(self, host, port, path):
        '''
        Constructor
        @param host: refers to the local host address
        @param port: port for the server to listen
        @param path: File system path where the data will be stored
        '''
        
        # server
        self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server.bind((host, port))
        self.server.listen(1)
        self.path = path
        self.app = None
    
    def setRXProperties(self, lo_off, fc, dec, gain, sync):
        '''
        Set the USRP RX properties
        
        @param lo_off: local oscillator offset (int)
        @param fc: Center frequency (float)
        @param dec: Decimation factor (int)
        @param gain: Gain of the receiver in dB (int)
        @param sync: True if the Hardware will use the GPSDO (boolean) 
        '''

        self.lo_off = lo_off
        self.fc = fc
        self.dec = dec
        self.gain = gain
        self.sync = sync
    
    def launch(self):
        '''
        calls startup
        '''
        print("i: launch SlaveRX")
        
        while True:
        
            sc, addr = self.server.accept()
            sc.settimeout(10*60)
            print("\n\ni: SlaveRX Connection from " + str(addr) + ", time " + time.strftime("%d-%m-%y/%H:%M:%S"))
            
            tic = time.time()
            
            try:
                self.__startup__(sc, addr)
            
            except Exception, e:
                print("e: " + str(e))
            
            sc.close()
            print("i: SlaveRX Connection closed, duration: " + str(time.time() - tic) + " [seg]\n\n")
            
        print("i: SlaveRX end launch")
        
    def record(self, prefix, at, signame):
        """
        @param prefix: prefix path folder where the signals are stored, e.g. /home/folder/
        @param at: attenuation factor  
        @param signame: filename of the signal 
        
        Start recording
        """
        # creating the folder
        folder = self.path + prefix
        folder = folder if (folder.endswith("/")) else folder + "/"
        Utils.ensure_dir(folder)
        
        # signal file
        filename = folder + signame + "_at" + str(at) +"_G" + str(self.gain) + ".dat"
        
        print("i: record filename = " + filename)
        
        self.app = RXApp(self.fc, self.dec, self.gain, "addr=192.168.10.2", self.sync, filename, self.lo_off)
        self.app.launch()
    
    def __startup__(self, sc, addr):
        '''
        Responsible for starting the application; for creating and showing
        the initial GUI.
        '''
        print("i: startup")
        
        msg = sc.recv(1024)
        if msg == "start":
            sc.send("ok")
            print("i: start ok")
            
            msg = sc.recv(1024)
            print("i: msg = " + msg)

            while msg != "finish":
                
                tic = time.time()
                
                if msg.find("startRec") >= 0:
                    # message "startRec:/prefix_path/:at:signame:"
                    print("i: startRec received")
                    values = msg.split(":")
                    
                    prefix = values[1]
                    at = float(values[2])
                    signame = values[3]
                    
                    self.record(prefix, at, signame)
                    sc.send("ok")
                
                elif msg.find("stopRec") >= 0:
                    print("i: stopRec received")
                    
                    if self.app.stopApp():
                        print("i: stopRec successful")
                        sc.send("ok")
                    else:
                        print("i: stopRec failed")
                        sc.send("error")
                        
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
        print "somebody call me!"
        self.__exit__()
        
        
if __name__ == '__main__':
    '''
    Creates an instance of the specified {@code Application}
    subclass, sets the {@code ApplicationContext} {@code
    application} property, and then calls the new {@code
    Application's} {@code startup} method.  The {@code launch} method is
    typically called from the Application's {@code main}:
    '''
    
    # Reading the properties
    confFile = "confRX.txt"
    if(len(sys.argv) > 1):
        arg = sys.argv[1]
        confFile = arg if len(arg) > 0 else confFile 
    else:
        print("working with default config file path")
    
    properties = readProperties(confFile)
    
    print("Properties:")
    for p in properties:
        print("\t" + p + " : " + properties[p])
       
    path = properties["rxpath"] 
    path = path if (path.endswith("/")) else path+"/"
    sync = True if properties["sync"] == "True" else False
    
    app = SlaveRX(properties["rxip"],
                  int(properties["rxport"]),
                  path)
    
    app.setRXProperties(int(properties["lo_off"]), 
                        float(properties["fc"]), 
                        int(properties["dec"]), 
                        int(properties["gain"]), 
                        sync)
    
    app.launch()
    exit()