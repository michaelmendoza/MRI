# ------------ Tissue (ms) -----------
# Gray Matter  : T1 - 920, T2 - 100
# White Matter : T1 - 790, T2 - 92
# Fat:           T1 - 270, T2 - 85
# Muscle:        T1 - 870, T2 - 47
# ------------ Tissue (ms) ------------
# Fat Water Spin Separation = 228 Hz

import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.pylab as plt
import Image

class TissueElement:
    def __init__(self, T1=0, T2=0, protonDensity=0, f0=0, Name=''):
        self.T1 = T1
        self.T2 = T2
        self.protonDensity = protonDensity
        self.f0 = f0
        self.Name = Name
        if T1 == 0 and T2 == 0:
            self.Name = 'Empty'

    def __str__(self):
        return "%s Element: T1=%f, T2=%f, ProtonDensity=%f, f0=%f" % (self.Name, self.T1, self.T2,
                                                                      self.protonDensity, self.f0)

class TissueElements:
    def __init__(self):
        self.tissueTypes = ['Empty', 'WhiteMatter', 'Fat']

    def getTissueElement(self, type='Empty', protonDensity=1):
        if type not in self.tissueTypes:
            return None

        Empty = TissueElement()
        WhiteMatter = TissueElement(790.0/1000.0, 92.0/1000.0, protonDensity, 0.0, 'WhiteMatter')
        Fat = TissueElement(270.0/1000.0, 85.0/1000.0, protonDensity, -428.0, 'Fat')                # -428 Hz @ 3T

        tissueTypes = {'Empty': Empty, 'WhiteMatter': WhiteMatter, 'Fat': Fat}
        return tissueTypes[type]

class TissueOptions:
    def __init__(self,
                 useFoPSF=True,      f0SD=1,
                 useRGradient=False, maxGradFreq=200,
                 noiseSD=0.0):
        self.useFoPSF = useFoPSF                    # Use Off-Resonance Point Spread Function (PSF) to create Gaussian Fo Distributions
        self.useRGradient = useRGradient            # Use Linear Off-Resonance Gradient
        self.f0SD = f0SD                            # Off-Resonance Std. Dev. in Hz
        self.maxGradFreq = maxGradFreq              # Max Frequency for Off-Resonance Linear Gradient
        self.noiseSD = noiseSD                      # Signal Noise Std. Dev.
        pass

    def __str__(self):
        return "TissueOptions: \n UseFoPSF = %s, FoSD = %f \n UseGradient = %s, MaxGradFreq = %f \n NoiseSD = %f" % (self.useFoPSF, self.f0SD, self.useRGradient, self.maxGradFreq, self.noiseSD)

class Tissue:
    def __init__(self, fileIndex=0, tissueOptions=TissueOptions()):
        files = ['RectPhantom64Tissues.png', 'Phantom64Tissues.png', 'Phantom256Tissues.png']
        self.tissueTypes = {'Empty': 0, 'Tissue1': 255, 'Tissue2': 91 }

        self.tissueOptions = tissueOptions
        self.loadTissueMask(files[fileIndex])
        self.generateTissue()
        self.generateResonancePSF(self.tissueOptions)

    def loadTissueMask(self, filename):
        imageMask = Image.open(filename).convert("L")
        imageMask = np.asmatrix(imageMask)
        self.Nx = np.shape(imageMask)[0]
        self.Ny = np.shape(imageMask)[1]

        self.mask = np.empty_like(imageMask)
        for r in range(self.Ny):
            for c in range(self.Nx):
                self.mask[r, c] = imageMask[r, c]

    def showTissueMask(self):
        plt.imshow(self.mask)
        #plt.imshow(self.mask, cmap = cm.Greys_r)
        plt.show()

    def generateTissue(self):
        elements = TissueElements()
        self.tissue = [[elements.getTissueElement('Empty') for c in range(self.Nx)] for r in range(self.Ny)]
        for r in range(self.Ny):
            for c in range(self.Nx):
                if self.mask[r,c] == self.tissueTypes['Tissue1']:
                    self.tissue[r][c] = elements.getTissueElement('WhiteMatter')
                elif self.mask[r,c] == self.tissueTypes['Tissue2']:
                    self.tissue[r][c] = elements.getTissueElement('Fat')
                else:
                    pass

    def generateResonancePSF(self, tissueOptions):
        for r in range(self.Ny):
            for c in range(self.Nx):
                self.tissue[r][c].f0 = self.tissue[r][c].f0 + tissueOptions.f0SD * np.random.uniform()

    def __str__(self):
        return "Tissue: \n Nx = %f, Ny = %f" % (self.Nx, self.Ny)

class PulseParameters:
    def __init__(self, alpha=math.pi/3.0, phi=0, dphi=0, TR=10.0/1000.0, TE=5.0/1000.0, Nr=200):
        self.alpha = alpha                        # Tip Angle
        self.phi = phi                            # Phi
        self.dphi = dphi                          # dPhi
        self.TR = TR                              # TR - Time between Flips
        self.TE = TE                              # TE - Time from Flip to Sample
        self.Nr = Nr                              # Number of Repetitions

    def __str__(self):
        return "Pulse Parameters: alpha=%f, phi=%f, dphi=%f \n TR=%f, TE=%f Nr=%i" % (self.alpha, self.phi, self.dphi, self.TR, self.TE, self.Nr)

class SamplingParameters:
    pass


class GradientParameters:
    def __init__(self, FOVx=0.256,     FOVy=0.256,
                       dt=1e-3,        tau=1e-3,
                       B0=3.0,         gamma=3 * 42.58e6):
        self.FOVx = FOVx                          # Field of View (FOV) in meters
        self.FOVy = FOVy                          # Field of View (FOV) in meters
        self.dt = dt                              # dt - Sample Period
        self.tau = tau                            # tau - Time for Phase Encode
        self.B0 = B0                              # B0 Field (T)
        self.gamma = gamma                        # Gyromagentic Ratio (Hz/T)
        #self.generateParameters()                 # x, y, dx, dy, dkx, dky, kxMax, kyMax, GxMax, Gyi, GyMax


    def __str__(self):
        return "Pulse Parameters: \n FOXx=%f,  FOXy=%f \n dt=%f,    tau=%f \n B0=%f,    gamma=%f" % (self.FOVx, self.FOVy, self.dt, self.tau, self.B0, self.gamma)
