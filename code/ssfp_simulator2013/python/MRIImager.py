
import kivy
kivy.require('1.0.6')

from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.relativelayout import RelativeLayout
from kivy.uix.image import Image
from kivy.properties import ObjectProperty, StringProperty, OptionProperty, NumericProperty,ReferenceListProperty, BoundedNumericProperty
from kivy.config import Config

from math import pi
from MRIpy import MRISimulator
import Tissue
import WaterFatSeparator

class Imager(FloatLayout):
    ''' Imager class - Main GUI class '''
    
    # Widget Varibles
    mainImageDisplay = ObjectProperty(None)
    mainPanel = ObjectProperty(None)
    algoModeSelector = ObjectProperty(None)
    statusLabel = StringProperty('Status')
    
    # Properties
    algoMode = OptionProperty('None', options=('None', 'SSFP', 'WaterFat'))
    
    # Init
    def __init__(self, **kwargs):
        super(Imager, self).__init__(**kwargs)
        # Init GUI State Machine
        self.algoMode = 'None'


    # GUI State Handlers
    def selectAlgoMode(self):
        mode = self.algoModeSelector.text
        print 'Changing Algorithm Mode:', mode
        if mode == 'None':
            self.algoMode = 'None'
            self.mainPanel.clear_widgets()
        elif mode == 'SSFP':
            self.algoMode = 'SSFP'
            self.mainPanel.clear_widgets()
            panel = SSPFOptionsPanel()
            self.mainPanel.add_widget(panel)

            display = ImageDisplay()
            display.opacity = 0
            self.mainImageDisplay.add_widget(display)
            panel.image = display
        elif mode == 'Water Fat Separation':
            self.algoMode = 'WaterFat'
            self.mainPanel.clear_widgets()
            panel = WaterFatOptionsPanel()
            self.mainPanel.add_widget(panel)

            display = ImageDisplay()
            display.opacity = 0
            self.mainImageDisplay.add_widget(display)
            panel.image = display

        return

class ImageDisplay(kivy.uix.image.Image):
    pass

class SSPFOptionsPanel(RelativeLayout):

    image = ObjectProperty(None)

    def plotOffResonance(self):

        T1 = 790.0/1000.0; T2 = 92.0/1000.0;
        sim = MRISimulator()
        f0, Mc = sim.SSFP_OffResonanceSpectrum( M0=1,
                                                alpha=self.alpha.value,
                                                phi=self.phi.value,
                                                dphi=self.dphi.value,
                                                Nr=int(self.nr.value),
                                                TR=self.tr.value/1000.0,
                                                TE=self.te.value/1000.0,
                                                T1=790.0/1000.0,
                                                T2=92.0/1000.0,
                                                Ns=200,
                                                BetaMax=2 * pi,
                                                f0=0.0)
        sim.saveMagPhasePlot(f0, Mc, 'temp.png')
        self.image.source = 'temp.png'
        self.image.reload()
        self.image.opacity = 1

    def showTissue(self):
        pass

    def showSSPFTissueSequence(self):
        options = Tissue.TissueOptions()
        tissue = Tissue.Tissue(1, options)
        p = Tissue.PulseParameters()
        sim = MRISimulator()
        MRISimulator.useOptimized = True
        img = sim.SSFP_MxyImage(tissue, p)
        sim.saveMagPhaseImage(img, 'temp.png')
        #sim.saveMagPhaseImage2(img, img, 'temp.png')
        #sim.saveMagPhaseImage3(img, img, img, 'temp.png')

        self.image.source = 'temp.png'
        self.image.reload()
        self.image.opacity = 1

class WaterFatOptionsPanel(RelativeLayout):

    image = ObjectProperty(None)

    def plotOffResonance(self):
        separator = WaterFatSeparator.Separator(alpha=self.alpha.value,
                                                phi=self.phi.value,
                                                dphi=self.dphi.value,
                                                Nr=int(self.nr.value),
                                                TR=self.tr.value/1000.0,
                                                TE=self.te.value/1000.0)
        separator.saveSimulation('temp.png')

        self.image.source = 'temp.png'
        self.image.reload()
        self.image.opacity = 1

    def showTissue(self):
        pass

    def showDixonTissueSimulation(self):
        pass

    def showSSPFTissueSimulation(self):
        separator = WaterFatSeparator.Separator(alpha=self.alpha.value,
                                                phi=self.phi.value,
                                                dphi=self.dphi.value,
                                                Nr=int(self.nr.value),
                                                TR=self.tr.value/1000.0,
                                                TE=self.te.value/1000.0)
        separator.saveSimulation('temp.png')

        self.image.source = 'temp.png'
        self.image.reload()
        self.image.opacity = 1

class PropertySlider(BoxLayout):
    text = StringProperty('')
    min = NumericProperty(0.)
    max = NumericProperty(100.)
    range = ReferenceListProperty(min, max)
    step = BoundedNumericProperty(0, min=0)
    value = NumericProperty(0.)

    def valueChanged(self):
        self.value = self.theSlider.value

class ImagerApp(App):
    def build(self):
        return Imager()

if __name__ == '__main__':
    Config.set('graphics', 'width', 1200)
    Config.set('graphics', 'height', 600)
    ImagerApp().run()