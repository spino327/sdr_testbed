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

from gnuradio.wxgui import fftsink2, scopesink2
from gnuradio import gr
from gnuradio import blks2
from grc_gnuradio import wxgui
from gnuradio.wxgui import forms

import wx

class RXGui(gr.hier_block2):
    '''
    
    This class have several input ports, the first one is for the Antenna signal and the second is for
    the demodulator output
    
    This class construct the GUI for the application, flow graph 
    
         2             2              2
    --->(throttle)--->(resampler)--->(fft)
                            |    2
                            --->(scope)
                            
        2 fft = 1 for raw antenna samples and 1 for matched filter
        2 scope = 1 for raw antenna samples and 1 for matched filter

    '''

    def __init__(self, app, gain, fc, samp_rate, inter, dec):
        '''
        in:
            - app = object of type RXApp
            - gain = gain
            - samp_rate = sample rate in Hertz
            - inter = interpolation factor
            - dec = decimation factor
        '''
        gr.hier_block2.__init__(self, "RXGui",
                                gr.io_signature(2, 2, gr.sizeof_gr_complex),
                                gr.io_signature(0, 0, 0))
          
        # instance variables
        self.app = app
        self.gui = wxgui.top_block_gui("BPSK RX") 
        self.nb = self.__createNoteBook()
        
        # controls
        self.gainTextBox = forms.text_box(
                parent=self.gui.GetWin(),
                value=gain,
                callback=self.setGain,
                label="gain",
                converter=forms.float_converter(),
        )
        self.fcTextBox = forms.text_box(
                parent=self.gui.GetWin(),
                value=fc,
                callback=self.setFc,
                label="fc",
                converter=forms.float_converter(),
        )
        self.startButton = wx.Button(self.gui.GetWin(), label="Record")
        self.startButton.Bind(wx.EVT_BUTTON, self.startRecording)
        
        # adding the visual components to the notebook
        self.gui.Add(self.gainTextBox)
        self.gui.Add(self.fcTextBox)
        self.gui.Add(self.startButton)
        
        #EO Controls
        
        # for port 1 Antenna samples (COMPLEX)
        self.throttleAn = gr.throttle(gr.sizeof_gr_complex, samp_rate)
        
        # resampler Antenna
        if inter == 1 and dec == 1:
            self.resamplerAn = gr.multiply_const_vcc((1,))
            print("i: resamplerAn not need")
        else:
            self.resamplerAn = blks2.rational_resampler_ccc(
                interpolation=inter,
                decimation=dec,
                taps=None,
                fractional_bw=None,
            )
        
#        self.cmp2arg1 = gr.complex_to_arg()
        
        self.fftAn = fftsink2.fft_sink_c(
                self.nb.GetPage(0).GetWin(),
                baseband_freq=0,
                y_per_div=5,
                y_divs=10,
                ref_level=-40,
                sample_rate= inter*samp_rate/dec,
                fft_size=512,
                fft_rate=10,
                average=True,
                avg_alpha=0.1,
                title="FFT Plot Antenna",
                peak_hold=False,
        )
        self.scope_IQAn = scopesink2.scope_sink_c(
                self.nb.GetPage(1).GetWin(),
                title="Scope IQ Antenna",
                sample_rate = inter*samp_rate/dec,
                v_scale=0.001,
                t_scale=0.001,
                ac_couple=False,
                xy_mode=False,
                num_inputs=1,
        )
        # adding the visual components to the notebook
        self.nb.GetPage(0).Add(self.fftAn.win)
        self.nb.GetPage(1).Add(self.scope_IQAn.win)
        
        # for port 2 Matched filter (COMPLEX)
        self.throttleMF = gr.throttle(gr.sizeof_gr_complex, samp_rate)
        
        # resampler MF
        if inter == 1 and dec == 1:
            self.resamplerMF = gr.multiply_const_vcc((1,))
            print("i: resamplerMF not need")
        else:
            self.resamplerMF = blks2.rational_resampler_ccc(
                interpolation=inter,
                decimation=dec,
                taps=None,
                fractional_bw=None,
            )
        
#        self.cmp2arg1 = gr.complex_to_arg()
        
        self.fftMF = fftsink2.fft_sink_c(
                self.nb.GetPage(2).GetWin(),
                baseband_freq=0,
                y_per_div=5,
                y_divs=10,
                ref_level=-40,
                sample_rate= inter*samp_rate/dec,
                fft_size=512,
                fft_rate=10,
                average=True,
                avg_alpha=0.1,
                title="FFT Plot MF",
                peak_hold=False,
        )
        self.scope_IQMF = scopesink2.scope_sink_c(
                self.nb.GetPage(3).GetWin(),
                title="Scope IQ MF",
                sample_rate = inter*samp_rate/dec,
                v_scale=0.0005,
                t_scale=0.001,
                ac_couple=False,
                xy_mode=False,
                num_inputs=1,
        )
        # adding the visual components to the notebook
        self.nb.GetPage(2).Add(self.fftMF.win)
        self.nb.GetPage(3).Add(self.scope_IQMF.win)
        
        # end of MF
        
        self.__makeConnections()
        
    def __createNoteBook(self):
        '''
        creates the NoteBook
        '''
        n1 = wx.Notebook(self.gui.GetWin(), style=wx.NB_RIGHT)
        n1.AddPage(wxgui.Panel(n1), "fft Ant")
        n1.AddPage(wxgui.Panel(n1), "scopeIQ Ant")
        
        n1.AddPage(wxgui.Panel(n1), "fft MF")
        n1.AddPage(wxgui.Panel(n1), "scopeIQ MF", True)
        self.gui.Add(n1)
        
        return n1
        
    def __makeConnections(self):
        '''
        uses the method connect(src, des)
        '''
        
        #Port 1
        self.connect((self, 0), (self.throttleAn, 0))
        self.connect((self.throttleAn, 0), (self.resamplerAn, 0))
        self.connect((self.resamplerAn, 0), (self.fftAn, 0))
        self.connect((self.resamplerAn, 0), (self.scope_IQAn, 0))
#        self.connect((self.resamplerAn, 0), (self.cmp2arg1, 0))
#        self.connect((self.cmp2arg1, 0), (self.fftAn, 0))
#        self.connect((self.cmp2arg1, 0), (self.scope_IQAn, 0))
#        null_sink = gr.null_sink(gr.sizeof_gr_complex*1)
#        self.connect((self, 0), null_sink)
        
        #Port 2
        self.connect((self, 1), (self.throttleMF, 0))
        self.connect((self.throttleMF, 0), (self.resamplerMF, 0))
        self.connect((self.resamplerMF, 0), (self.fftMF, 0))
        self.connect((self.resamplerMF, 0), (self.scope_IQMF, 0))
#        self.connect((self.resamplerDem, 0), (self.cmp2arg2, 0))
#        self.connect((self.cmp2arg2, 0), (self.fftDem, 0))
#        self.connect((self.cmp2arg2, 0), (self.scope_IQDem, 0))
        
    def Run(self):
        '''
        calls the Run method in the gui object
        '''
        self.gui.Run(True)
    
    def setFc(self, fc):
        self.fcTextBox.set_value(fc)
        self.app.setFc(fc)
    
    def setGain(self, gain):
        self.gainTextBox.set_value(gain)
        self.app.setGain(gain)
        
    def startRecording(self, event):
        self.app.startRecording()
   
#if __name__ == "__main__":
#    
#    tb = gr.top_block()
#    signalRaw = gr.sig_source_c(1e4, gr.GR_SIN_WAVE, 350, 1)
#    signalDem = gr.sig_source_c(1e4, gr.GR_TRI_WAVE, 200, 1)
#    signalCL = gr.sig_source_c(1e4, gr.GR_SIN_WAVE, 350, 1)
#    signalAGC = gr.sig_source_c(1e4, gr.GR_TRI_WAVE, 200, 1)
#    temp = RXGui(None, 1, 0, 1e4, 1, 1)
#    tb.connect(signalRaw, (temp, 0))
#    tb.connect(signalAGC, (temp, 1))
#    tb.connect(signalCL, (temp, 2))
#    tb.connect(signalDem, (temp, 3))
#    tb.start()
#    temp.Run()
        
        