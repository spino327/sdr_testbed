'''
Copyright (c) 2010, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Dec, 2010
'''

#from gnuradio import usrp2
from gnuradio import uhd

import time

def getUhdUSRPSink(fc, lo_off, inter, gain, addr, sync):
    """
    Tx
    def getUhdUSRPSink(fc, inter, gain, addr, sync)
    
    in: 
    - fc = center frequency
    - lo_off = LO off
    - inter = interpolation factor
    - gain = gain in the tx, only with 2450
    - addr = ip address, format = "addr=ip, mimo_mode="
    - sync = if True them sync with external clock
    
    out:
    (usrp2, baseband_freq, dxc_freq)
    - usrp2 sink object
    - baseband_freq
    - dxc_freq
    """
    
    rf_sink = uhd.usrp_sink(addr, uhd.io_type_t.COMPLEX_FLOAT32, 1)
    
    # gain
    gRange = rf_sink.get_gain_range()
    
    if gRange.start() <> gRange.stop():
        rf_sink.set_gain(gain)
        print("i: set_gain = " + str(rf_sink.set_gain()))
    else:
        print("i: this daughterboard not support the set gain behavior")
    
    # samo freq
    rf_sink.set_samp_rate(rf_sink.get_clock_rate()/inter)
    print("i: samp_freq = " + str(rf_sink.get_samp_rate()))
    
    # center freq
    freqRange = rf_sink.get_freq_range()
    
    if float(freqRange.start()) > fc and float(freqRange.stop()) < fc:
        fc = float(freqRange.start()+freqRange.stop())/2
        print("e: fc have to be between [" + str(freqRange.start()) + ", " + str(freqRange.stop()) + "]")
        
        
    tune_request = uhd.tune_request_t(fc, lo_off)
    tune_res = rf_sink.set_center_freq(tune_request)
#    tune_res = rf_sink.set_center_freq(fc)
        
    rf_freq = tune_res.actual_rf_freq
    dxc_freq = 0
    print("i: actual_rf_freq = " + str(rf_freq) + 
          "\ni: actual_dsp_freq = " + str(tune_res.actual_dsp_freq) + 
          "\ni: target_rf_freq = " + str(tune_res.target_rf_freq) +
          "\ni: target_dsp_freq = " + str(tune_res.target_dsp_freq))
    
    # sync
    if type(sync) == type(False) and bool(sync):
#        print "i: sync to pps = " + str(rf_sink.set_time_unknown_pps(uhd.time_spec_t()))
#        print "i: sync to pps = " + str(rf_sink.set_time_next_pps(uhd.time_spec_t()))
       
        # Common references signals
        ccfg = uhd.clock_config_t()
        ccfg.ref_source = uhd.clock_config_t.REF_SMA
        ccfg.pps_source = uhd.clock_config_t.PPS_SMA
        ccfg.pps_polarity = uhd.clock_config_t.PPS_NEG
        
        rf_sink.set_clock_config(ccfg)
        #EO CRS
        
        # sync the device time
#        last_pps_time = rf_sink.get_time_last_pps()
#        while (last_pps_time != rf_sink.get_time_last_pps()):
#            time.sleep(0.1)
            
        rf_sink.set_time_next_pps(uhd.time_spec_t(0.0))
        print("i: sync == True")
        
    else:
        print("i: sync == False")
        
#    time.sleep(2)
    
    return (rf_sink, rf_freq, dxc_freq)

def getUhdUSRP2Source(fc, lo_off, dec, gain, addr, sync):
    """
    Rx
    def getUhdUSRP2Source(fc, dec, gain, addr, sync):
    
    in: 
    - fc = center frequency
    - dec = decimation factor
    - gain = gain in the tx, only with 2450
    - addr = ip address, format = "addr=ip, mimo_mode="
    - sync = if True them sync with external clock
    
    out:
    (usrp2, baseband_freq, dxc_freq)
    - usrp2 source object
    - baseband_freq
    - dxc_freq
    """
    
    rf_src = uhd.usrp_source(addr, uhd.io_type_t.COMPLEX_FLOAT32, 1)
    
    # gain
    gRange = rf_src.get_gain_range()
    
    if gRange.start() <> gRange.stop():
        rf_src.set_gain(gain)
        print("i: set_gain = " + str(rf_src.get_gain()))
    else:
        print("i: this daughterboard not support the set gain behavior")
    
    # samo freq
    rf_src.set_samp_rate(rf_src.get_clock_rate()/dec)
    print("i: samp_freq = " + str(rf_src.get_samp_rate()))
    
    # center freq
    freqRange = rf_src.get_freq_range()
    
    if float(freqRange.start()) > fc and float(freqRange.stop()) < fc:
        fc = float(freqRange.start()+freqRange.stop())/2
        print("e: fc have to be between [" + str(freqRange.start()) + ", " + str(freqRange.stop()) + "]")
    
    tune_request = uhd.tune_request_t(fc, lo_off)
    tune_res = rf_src.set_center_freq(tune_request)
        
    rf_freq = tune_res.actual_rf_freq
    dxc_freq = 0
    print("i: actual_rf_freq = " + str(rf_freq) + 
          "\ni: actual_dsp_freq = " + str(tune_res.actual_dsp_freq) + 
          "\ni: target_rf_freq = " + str(tune_res.target_rf_freq) +
          "\ni: target_dsp_freq = " + str(tune_res.target_dsp_freq))
    
    # sync
    if type(sync) == type(False) and bool(sync):
##        print "i: sync to pps = " + str(rf_src.set_time_unknown_pps(uhd.time_spec_t()))
##        print "i: sync to pps = " + str(rf_src.set_time_next_pps(uhd.time_spec_t()))
#        print "i: sync to pps = " + str(rf_src.set)
#        rf_src.set_time_next_pps(uhd.time_spec_t(uhd.time_spec_t.get_system_time().get_real_secs()+1))
        
        # Common references signals
        ccfg = uhd.clock_config_t()
        ccfg.ref_source = uhd.clock_config_t.REF_SMA
        ccfg.pps_source = uhd.clock_config_t.PPS_SMA
        ccfg.pps_polarity = uhd.clock_config_t.PPS_NEG
        
        rf_src.set_clock_config(ccfg)
        #EO CRS
        
        # sync the device time
#        last_pps_time = rf_src.get_time_last_pps()
#        while (last_pps_time != rf_src.get_time_last_pps()):
#            time.sleep(0.1)
#        rf_src.set_time_unknown_pps(uhd.time_spec_t(0.0))
        
        rf_src.set_time_next_pps(uhd.time_spec_t(0.0))
        print("i: sync == True")
        
        # EO SDT
    else:
        print("i: sync == False")   
        
#    time.sleep(1)
        
    return (rf_src, rf_freq, dxc_freq)
    
    
#def getUSRP2Sink(fc, inter, gain, eth, sync):
#    """
#    Tx
#    def getUSRP2Sink(fc, inter, gain, eth, sync):
#    
#    in: 
#    - fc = center frequency
#    - inter = interpolation factor
#    - gain = gain in the tx, only with 2450
#    - eth = ethernet interface name(String)
#    - sync = if True them sync with external clock
#    
#    out:
#    (usrp2, baseband_freq, dxc_freq)
#    - usrp2 sink object
#    - baseband_freq
#    - dxc_freq
#    """
#    rf_sink = usrp2.sink_32fc(eth)
#    
#    if rf_sink.gain_max() <> rf_sink.gain_min():
#        print("i: set_gain = " + str(rf_sink.set_gain(gain)))
#    else:
#        print("i: this daughterboard not support the set gain behavior")
#    
#    print("i: set inter = " +str(rf_sink.set_interp(inter)))
#    tune_res = rf_sink.set_center_freq(fc)
#    
#    baseband_freq = tune_res.baseband_freq
#    dxc_freq = tune_res.dxc_freq
#    print("i: set center freq = " + str(baseband_freq) + 
#          "\ni: residual freq = " + str(tune_res.residual_freq) + 
#          "\ni: dxc freq = " + str(dxc_freq))
#    
#    if type(sync) == type(False) and bool(sync):
##        print "i: sync to mimo = " + str(rf_sink.config_mimo(usrp2.MC_WE_LOCK_TO_MIMO))
#        print "i: sync to pps = " + str(rf_sink.sync_to_pps())
#        print "i: sync to ref clock = " + str(rf_sink.config_mimo(usrp2.MC_WE_LOCK_TO_SMA))
#    else:
#        print "i: ref clock don't lock = " + str(rf_sink.config_mimo(usrp2.MC_WE_DONT_LOCK))
#    
#    print("i: samp_freq = " + str(float(rf_sink.dac_rate())/inter))
#    
#    return (rf_sink, baseband_freq, dxc_freq)

#def getUSRP2Source(fc, dec, gain, eth, sync):
#    """
#    Rx
#    def getUSRP2Source(fc, dec, gain, eth, sync):
#    
#    in: 
#    - fc = center frequency
#    - dec = decimation factor
#    - gain = gain in the tx, only with 2450
#    - eth = ethernet interface name(String)
#    - sync = if True them sync with external clock
#    
#    out:
#    (usrp2, baseband_freq, dxc_freq)
#    - usrp2 source object
#    - baseband_freq
#    - dxc_freq
#    """
#    
#    rf_src = usrp2.source_32fc(eth)
#    
#    print("i: set_gain = " + str(rf_src.set_gain(gain)))
#    print("i: set int = " +str(rf_src.set_decim(dec)))
#    tune_res = rf_src.set_center_freq(fc)
#    
#    baseband_freq = tune_res.baseband_freq
#    dxc_freq = tune_res.dxc_freq
#    print("i: set center freq = " + str(baseband_freq) + 
#          "\ni: residual freq = " + str(tune_res.residual_freq) + 
#          "\ni: dxc freq = " + str(dxc_freq))
#    
#    if type(sync) == type(False) and bool(sync):
##        print "i: sync to mimo = " + str(rf_src.config_mimo(usrp2.MC_PROVIDE_CLK_TO_MIMO))
#        print "i: sync to pps = " + str(rf_src.sync_to_pps())
#        print "i: sync to ref clock = " + str(rf_src.config_mimo(usrp2.MC_WE_LOCK_TO_SMA))
#    else:
#        print "i: ref clock don't lock = " + str(rf_src.config_mimo(usrp2.MC_WE_DONT_LOCK))
#    
#    print("i: samp_freq = " + str(float(rf_src.adc_rate())/dec))
#    
#    return (rf_src, baseband_freq, dxc_freq)
        