import numpy as np

import matplotlib.cm as cm
from matplotlib.pylab import plt

import MRIpy
import Tissue
import math

class Separator():

    # Initalization
    def __init__(self, alpha, phi, dphi, TR, TE, Nr, mode='SSFP'):
        '''Create a Separator. Sets alpha, phi, dphi, TR, TE, Nr parameters. Loads a tissue model.'''
        self.p = None
        self.mode = mode
        self.simulation = None

        self.loadTissue()
        self.generateParameters(alpha, phi, dphi, TR, TE, Nr)

    def loadTissue(self):
        '''Loads a tissue model to Separator, and loads tissue parameters to Separator'''
        # Load Tissue
        options = Tissue.TissueOptions()
        self.tissue = Tissue.Tissue(1, options)

        # Load Tissues in use
        self.tissues = []
        elements = Tissue.TissueElements()
        self.tissues.append(elements.getTissueElement('WhiteMatter'))
        self.tissues.append(elements.getTissueElement('Fat'))

    def generateParameters(self, alpha, phi, dphi, TR, TE, Nr):
        '''Generates a set of pulse parameters for Water/Fat Separation. For SSFP mode this generates 4 sets of parameters at different
            TE and dPhi'''
        elements = Tissue.TissueElements()
        df = abs(elements.getTissueElement('WhiteMatter').f0 - elements.getTissueElement('Fat').f0)

        p0 = Tissue.PulseParameters(alpha, phi, dphi, TR, TE - 2 / (2 * df), Nr)
        p1 = Tissue.PulseParameters(alpha, phi, dphi, TR, TE - 1 / (2 * df), Nr)
        p2 = Tissue.PulseParameters(alpha, phi, dphi, TR, TE,                Nr)
        p3 = Tissue.PulseParameters(alpha, phi, dphi, TR, TE + 1 / (2 * df), Nr)
        p4 = Tissue.PulseParameters(alpha, phi, dphi, TR, TE + 2 / (2 * df), Nr)
        self.p = [p0, p1, p2, p3, p4]

    # Tissue Spectra
    def saveTissueSpectra(self, filename):
        if self.p is None:
            return

        sim = MRIpy.MRISimulator()
        for n in range(len(self.tissues)):
            p = self.p[0]
            t = self.tissues[n]
            f0, Mc = sim.SSFP_OffResonanceSpectrum(M0=1, alpha=p.alpha, phi=p.phi, dphi=p.dphi, Nr=p.Nr, TR=p.TR,
                                                   TE=p.TE, T1=t.T1, T2=t.T2, Ns=200, BetaMax=2*math.pi, f0=t.f0)
            sim.saveMagPhasePlot(f0, Mc, filename)

    # Separation Simulation
    def separationSimulation(self):
        if self.mode is 'SSFP':
            pass
        elif self.mode is 'Dixon':
            sim = MRIpy.MRISimulator()
            img = sim.SSFP_MxyImage(self.tissue, self.p[0])
            img2 = sim.SSFP_MxyImage(self.tissue, self.p[1])
            img3 = sim.SSFP_MxyImage(self.tissue, self.p[2])

            self.simulation = Dixon()
            self.simulation.threePointDixon(img, img2, img3)

    def saveSimulation(self, filename):
        if self.mode is 'SSFP':
            pass
        elif self.mode is 'Dixon':
            self.simulation.saveWaterFatImage(filename)
        pass


class Dixon():

    def __init__(self):
        self.water = None
        self.fat = None

    def twoPointDixon(self, image1, image2):
        self.water = 0.5 * (image1 + image2)
        self.fat = 0.5 * (image1 - image2)
        return self.water, self.fat

    def threePointDixon(self, image1, image2, image3):
        phi = np.angle(np.conj(image1) * image3) / 2.0;
        self.water = 0.5 * (image1 + image2 * np.exp(-1j * phi))
        self.fat = 0.5 * (image1 - image2 * np.exp(-1j * phi))
        return self.water, self.fat

    def saveWaterFatImage(self, filename):
        mag = np.absolute(self.water).astype('float')
        phase = np.angle(self.water).astype('float')
        mag2 = np.absolute(self.fat).astype('float')
        phase2 = np.angle(self.fat).astype('float')

        plt.subplot(221)
        plt.imshow(mag, cmap=cm.Greys_r)
        plt.title('Water Magnitude')
        plt.axis('off')
        plt.subplot(222)
        plt.imshow(phase, cmap=cm.Greys_r)
        plt.title('Water Phase')
        plt.axis('off')
        plt.subplot(223)
        plt.imshow(mag2, cmap=cm.Greys_r)
        plt.title('Fat Magnitude')
        plt.axis('off')
        plt.subplot(224)
        plt.imshow(phase2, cmap=cm.Greys_r)
        plt.title('Fat Phase')
        plt.axis('off')
        plt.show()

    def saveWaterFatImage(self, filename):
        mag = np.absolute(self.water).astype('float')
        phase = np.angle(self.water).astype('float')
        mag2 = np.absolute(self.fat).astype('float')
        phase2 = np.angle(self.fat).astype('float')

        plt.subplot(221)
        plt.imshow(mag, cmap=cm.Greys_r)
        plt.title('Water Magnitude')
        plt.axis('off')
        plt.subplot(222)
        plt.imshow(phase, cmap=cm.Greys_r)
        plt.title('Water Phase')
        plt.axis('off')
        plt.subplot(223)
        plt.imshow(mag2, cmap=cm.Greys_r)
        plt.title('Fat Magnitude')
        plt.axis('off')
        plt.subplot(224)
        plt.imshow(phase2, cmap=cm.Greys_r)
        plt.title('Fat Phase')
        plt.axis('off')
        plt.savefig(filename)

def DixonTest():
    pass

def SSPFSeparationTest():
    pass