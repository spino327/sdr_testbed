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
from grc_gnuradio.wxgui.top_block_gui import default_gui_size

import wx

class TXGui(gr.hier_block2):
    '''
    This class construct the GUI for the application, flow graph
      
      
    --->(throttle)--->(resampler)--->(FFT)
                            |
                            --->(scope)
                            |
                            |
                            --->(complex to arg)--->(scope)
                                        |
                                        |
                                        --->(FFT)
      
    '''

    def __init__(self, app, samp_rate, inter, dec):
        '''
        in:
            - app = object of type TXApp
            - samp_rate = sample rate in Herz
            - inter = interpolation factor
            - dec = decimation factor
        '''
        gr.hier_block2.__init__(self, "TXGui",
                                gr.io_signature(1, 1, gr.sizeof_gr_complex),
                                gr.io_signature(0, 0, 0))
          
        # instance variables
        self.app = app
        self.gui = wxgui.top_block_gui("BPSK TX", size=default_gui_size) 
        self.nb = self.__createNoteBook()
        
        self.throttle = gr.throttle(gr.sizeof_gr_complex, samp_rate)
        
        # resampler
        if inter == 1 and dec == 1:
            self.resampler = gr.multiply_const_vcc((1,))
            print("i: resampler not need")
        else:
            self.resampler = blks2.rational_resampler_ccc(
                interpolation=inter,
                decimation=dec,
                taps=None,
                fractional_bw=None,
            )
        # samples
        self.fft = fftsink2.fft_sink_c(
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
                title="FFT Plot",
                peak_hold=False,
        )
        self.scope_IQ = scopesink2.scope_sink_c(
                self.nb.GetPage(1).GetWin(),
                title="Scope IQ",
                sample_rate = inter*samp_rate/dec,
                v_scale=0.001,
                t_scale=0.001,
                ac_couple=False,
                xy_mode=False,
                num_inputs=1,
        )
        
        # adding the visual components to the notebook
        self.nb.GetPage(0).Add(self.fft.win)
        self.nb.GetPage(1).Add(self.scope_IQ.win)
        
        # phase
        self.cmp2arg = gr.complex_to_arg()
        
        self.fftPhase = fftsink2.fft_sink_f(
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
                title="FFT phase Plot",
                peak_hold=False,
        )
        self.scopePhase = scopesink2.scope_sink_f(
                self.nb.GetPage(3).GetWin(),
                title="Scope phase",
                sample_rate = inter*samp_rate/dec,
                v_scale=0.001,
                t_scale=0.001,
                ac_couple=False,
                xy_mode=False,
                num_inputs=1,
        )
        # adding the visual components to the notebook
        self.nb.GetPage(2).Add(self.fftPhase.win)
        self.nb.GetPage(3).Add(self.scopePhase.win)
        
        self.__makeConnections()
        
    def __createNoteBook(self):
        '''
        creates the NoteBook
        '''
        n1 = wx.Notebook(self.gui.GetWin(), style=wx.NB_RIGHT)
        n1.AddPage(wxgui.Panel(n1), "fft")
        n1.AddPage(wxgui.Panel(n1), "scope IQ")
        n1.AddPage(wxgui.Panel(n1), "fft phase")
        n1.AddPage(wxgui.Panel(n1), "scope phase")
        self.gui.Add(n1)
        
        return n1
        
    def __makeConnections(self):
        '''
        uses the method connect(src, des)
        '''
        self.connect(self, (self.throttle, 0))
        self.connect((self.throttle, 0), (self.resampler, 0))
        self.connect((self.resampler, 0), (self.fft, 0))
        self.connect((self.resampler, 0), (self.scope_IQ, 0))
        
        self.connect((self.resampler, 0), (self.cmp2arg, 0))
        self.connect((self.cmp2arg, 0), (self.fftPhase, 0))
        self.connect((self.cmp2arg, 0), (self.scopePhase, 0))
        
        
    def Run(self):
        '''
        calls the Run method in the gui object
        '''
        self.gui.Run(True)
        
    
   
#if __name__ == "__main__":
#    
#    tb = gr.top_block()
#    signal = gr.sig_source_c(1e3, gr.GR_SIN_WAVE, 350, 1)
#    temp = TXGui(None, 1e3, 1, 1)
#    tb.connect(signal, temp)
#    tb.start()
#    temp.Run()