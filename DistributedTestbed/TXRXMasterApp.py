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
from util.PropertyReader import readProperties

class TXRXMasterApp(object):
    '''
    Application life cycle
    '''

    def __init__(self, slaveTXIP, slaveTXPort, slaveRXIP, slaveRXPort, sigs, ats, timeSine, timeSig, txPrefix, rxPrefix):
        '''
        Connect to the instances of slaveTX and slaveRX:

        @param slaveTXIP : string
        @param slaveTXPort : int
        @param slaveRXIP : String
        @param slaveRXPort : int
        @param sigs : signals to transmit
        @param ats: attenuation to use with each signal
        @param timeSine : int
        @param timeSig : int
        @param txPrefix : prefix for the transmitted data at TX
        @param rxPrefix : prefix for the received data at RX
        '''
        
        if type(slaveTXIP) != type("1") or type(slaveRXIP) != type("1") or type(slaveTXPort) != type(1) or type(slaveRXPort) != type(1):
            raise Exception("e:Parameter error")
        
        # signals
        self.timeSine = timeSine
        self.timeSig = timeSig
        
        #creating the sokects
        self.tx = None
        self.slaveTXIP = slaveTXIP
        self.slaveTXPort = slaveTXPort
        self.txPrefix = txPrefix if (txPrefix.endswith("/")) else txPrefix + "/"
        
        self.rx = None
        self.slaveRXIP = slaveRXIP
        self.slaveRXPort = slaveRXPort
        self.rxPrefix = rxPrefix if (rxPrefix.endswith("/")) else rxPrefix + "/"
        
        # position in the vector of attenuations
        self.sigs = sigs if len(sigs) == len(ats) else sigs*len(ats)
        self.ats = ats
        
         
    def launch(self):
        '''
        calls startup
        '''
        print("i: launch")
        
        init = time.time()
#        for at in self.ats:
        for pos in range(len(ats)):
            
            time.sleep(2)

            print("\n\ni: Processing at: " + str(self.ats[pos]) + " and signal: " + str(self.sigs[pos]))
            tic = time.time()
            
            try:
                self.process(self.ats[pos], self.sigs[pos])
            
            except Exception, e:
                print("e: " + str(e))
                
            print("i: End processing at: " + str(self.ats[pos]) + " and signal: " + str(self.sigs[pos]) + "\n\n")
    
        print("i: End TXRXMasterApp, duration: " + str(time.time() - init) + " [seg]")
    
    def process(self, at, sig):
        '''
        @param at: atenuation factor
        '''
        
        socket_to = 10*60# seg
        self.tx = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.tx.settimeout(socket_to)
        
        self.rx = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.rx.settimeout(socket_to)
        
        print("i: TX socket time out set to " + str(self.tx.gettimeout()))
        self.tx.connect((self.slaveTXIP, self.slaveTXPort))
        print("i: TX connection successful")
        
        print("i: RX socket time out set to " + str(self.rx.gettimeout()))
        self.rx.connect((self.slaveRXIP, self.slaveRXPort))
        print("i: RX connection successful")
        
        # startTx - ok
        tic = time.time()
        self.tx.send("start")
        resTx = self.tx.recv(1024)
        print("i: TX start cmd duration: " + str(time.time() - tic) + " [seg]")
        
        # startRx - ok
        tic = time.time()
        self.rx.send("start")
        resRx = self.rx.recv(1024)
        print("i: RX start cmd duration: " + str(time.time() - tic) + " [seg]")
        
        if resTx == "ok" and resRx == "ok":
            
            # set Signal
            tic = time.time()
            if self.sendSetSignalCmd(sig):
                print("i: TX set signal cmd duration: " + str(time.time() - tic) + " [seg]")
                
                # startRec Sine
                tic = time.time()
                if self.sendStartRecordingCmd(self.rxPrefix + "fading/", 1, sig):
                    print("i: RX [Sine] start recording cmd duration: " + str(time.time() - tic) + " [seg]")
                
                    # Sending Sine
                    tic = time.time()
                    if self.sendSineCmd():
                        print("i: TX send sine cmd duration: " + str(time.time() - tic) + " [seg]")
                    else:
                        print("e: TX send sine couldn't be executed by SlaveTX")
                    
                    # StopRec Sine
                    time.sleep(1)
                    tic = time.time()
                    if self.sendStopRecordingCmd():
                        print("i: [Sine] stop recording cmd duration: " + str(time.time() - tic) + " [seg]")
                    else:
                        print("e: RX SlaveRX couldn't stop recording [sine]")
                        
                else:
                    print("e: RX start recording couldn't be executed by SlaveRX")
                print("")
                time.sleep(1)
                
#                # Carrier recovery - Frequency synchronization
#                tic = time.time()
#                if self.sendStartRecordingCmd(self.rxPrefix + "carrier/", 1, sig):
#                    print("i: RX [Carrier recovery] start recording cmd duration: " + str(time.time() - tic) + " [seg]")
#                
#                    # Sending Carrier recovery
#                    tic = time.time()
#                    if self.sendCarrierRecoveryCmd():
#                        print("i: TX send Carrier recovery cmd duration: " + str(time.time() - tic) + " [seg]")
#                    else:
#                        print("e: TX send Carrier recovery couldn't be executed by SlaveTX")
#                    
#                    # StopRec Carrier recovery
#                    time.sleep(1)
#                    tic = time.time()
#                    if self.sendStopRecordingCmd():
#                        print("i: [Carrier recovery] stop recording cmd duration: " + str(time.time() - tic) + " [seg]")
#                    else:
#                        print("e: RX SlaveRX couldn't stop recording [Carrier recovery]")
#                        
#                else:
#                    print("e: RX start recording couldn't be executed by SlaveRX")
#                print("")
#                time.sleep(1)
                
                # StartRec Signal
                tic = time.time()
                if self.sendStartRecordingCmd(self.rxPrefix, at, sig):
                    print("i: [Signal] start recording cmd duration: " + str(time.time() - tic) + " [seg]")
            
                    # Sending signal
                    tic = time.time()
                    
                    # delay for get more samples of the noise
                    time.sleep(0.5);
                    
                    if self.sendSignalCmd(at):
                        print("i: TX send Signal cmd duration: " + str(time.time() - tic) + " [seg]")
                    else:
                        print("e: TX send signal couldn't be executed by SlaveTX")
                        
                    # Stop Recording
                    time.sleep(2)
                    print("i: RX [Signal] sending stop recording cmd")
                    tic = time.time()
                    if self.sendStopRecordingCmd():
                        print("i: RX [Signal] stop recording cmd duration: " + str(time.time() - tic) + " [seg]")
                    else:
                        print("e: RX SlaveRX couldn't stop recording [Signal]")
                
                else:
                    print("e: RX SlaveRX couldn't start recording [Signal]")
                 
            else:
                print("e: TX SlaveTX couldn't set signal")
                
        else:
            print("e: TX and RX Slaves didn't response to start command")
                
        # finish anyway
        self.tx.send("finish")
        self.tx.close()
        print("i: TX end communication")

        self.rx.send("finish")
        self.rx.close()
        print("i: RX end communication")
    
    def sendSetSignalCmd(self, sig):
        """
        @param sig: signal file name 
        
        Send the set signal (path) command to SlaveTX
        
        message    "setSig:TXPrefix/signal_name:"
        """
        self.tx.send("setSig:" + self.txPrefix + sig + ":")
        response = self.tx.recv(1024)
        
        if response == "ok":
            return True
        
        return False
    
    def sendSineCmd(self):
        """
        Send the send Sine command to SlaveTX
        
        message    "sendSine:timeSine:"
        """
        self.tx.send("sendSine:" + str(self.timeSine) + ":")
        response = self.tx.recv(1024)
        
        if response == "ok":
            return True
        
        return False
    
    def sendCarrierRecoveryCmd(self):
        """
        Send the send Carrier Recovery command to SlaveTX
        
        message    "sendCR"
        """
        self.tx.send("sendCR")
        response = self.tx.recv(1024)
        
        if response == "ok":
            return True
        
        return False
    
    def sendSignalCmd(self, at):
        """
        @param at: attenuation factor to be used in the transmission
        
        Send the send Signal command to SlaveTX
        
        message:    "sendSignal:timeSig:at:"
        """
        self.tx.send("sendSignal:" + str(self.timeSig) + ":" + str(at) + ":")
        response = self.tx.recv(1024)
        
        if response == "ok":
            return True
        
        return False
    
    def sendStartRecordingCmd(self, prefix, at, signame):
        """
        @param prefix: prefix path folder where the signals are stored, e.g. /home/folder/
        @param at: attenuation factor  
        @param signame: filename of the signal 
        
        Send start recording command to SlaveRX
        
        message:    "startRec:/prefix_path/:at:signame:"
        """
        # fixing the signame by removing the file extension (if any)
        signame = signame[:signame.rfind(".")]
        
        self.rx.send("startRec:" + prefix + ":" + str(at) + ":" + signame + ":")
        response = self.rx.recv(1024)
        
        if response == "ok":
            return True
        
        return False
    
    def sendStopRecordingCmd(self):
        """
        Send stop recording command to SlaveRX
        
        message:    "stopRec"
        """
        self.rx.send("stopRec")
        response = self.rx.recv(1024)
        
        if response == "ok":
            return True
        
        return False
    
    def __startup__(self):
        '''
        Responsible for starting the application; for creating and showing
        the initial GUI.
        '''
        print("i: startup")
        
        self.slaveTX.start()
        time.sleep(1)
        self.slaveRX.start()
        
        self.slaveTX.join()
        self.slaveRX.join()
        
        print("i: end startup")
        self.__exit__()
       
        
    def __exit__(self):
        '''
        This method runs on the event dispatching thread.
        '''
        print("somebody call me!")
        

if __name__ == '__main__':
    '''
    Creates an instance of the specified {@code Application}
    subclass, sets the {@code ApplicationContext} {@code
    application} property, and then calls the new {@code
    Application's} {@code startup} method.  The {@code launch} method is
    typically called from the Application's {@code main}:
    '''
    
    # Reading the properties
    confFile = "confMaster.txt"
    if(len(sys.argv) > 1):
        arg = sys.argv[1]
        confFile = arg if len(arg) > 0 else confFile 
    else:
        print("working with default config file path")
    
    properties = readProperties(confFile)
    
    print("Properties:")
    for p in properties:
        print("\t" + p + " : " + properties[p])

    # attenuations vector
    input = properties["txat"]
    if len(input) > 0:
        ats = input.split(";")
        for i in range(len(ats)):
            ats[i] = float(ats[i].strip())
            
    # signals vector
    input = properties["txsig"]
    if len(input) > 0:
        sigs = input.split(";")
        
        # signals number should be either 1 or the length of the ats vector
        if len(sigs) != 1 and len(sigs) != len(ats):
            raise Exception("e:Number of signals and attenuation error")

    # Creating the app    
    app = TXRXMasterApp(properties["txip"], int(properties["txport"]), 
                        properties["rxip"], int(properties["rxport"]), 
                        sigs, 
                        ats, 
                        int(properties["tsine"]), 
                        int(properties["tsig"]),
                        properties["txprefix"], 
                        properties["rxprefix"])
    
    app.launch()
    exit()