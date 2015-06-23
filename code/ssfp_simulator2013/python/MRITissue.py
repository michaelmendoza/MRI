from math import *
from cmath import *
from numpy import *
import matplotlib.cm as cm
from matplotlib.pylab import plt
from mpl_toolkits.mplot3d import Axes3D
import Image

class TissueElement:
    def __init__(self, T1 = 0, T2 = 0, protonDensity = 0, f0 = 0, Name = ''):
        self.T1 = T1
        self.T2 = T2
        self.protonDensity = protonDensity
        self.f0 = f0
        self.Name = Name
        if T1 == 0 and T2 == 0:
            self.Name = 'Empty'
    
    def __str__(self):
        return "%s Element: T1=%f, T2=%f, ProtonDensity=%f, f0=%f" % (self.Name, self.T1, self.T2, self.protonDensity, self.f0)

class TissueElements:
    def __init__(self):
        self.tissueTypes = ['Empty', 'WhiteMatter', 'Fat']

    def getTissueElement(self, type = 'Emtpy', protonDensity = 1):
        if type not in self.tissueTypes:
            return None
        
        Empty = TissueElement()
        WhiteMatter = TissueElement(790.0/1000.0, 92.0/1000.0, protonDensity, 0.0, 'WhiteMatter')
        Fat = TissueElement(270.0/1000.0, 85.0/1000.0, protonDensity, -428.0, 'Fat')                # -428 Hz @ 3T
        
        tissueTypes = {'Empty': Empty, 'WhiteMatter': WhiteMatter, 'Fat': Fat}gt
        return tissueTypes[type]

class TissueOptions:
    def __init__(self,
                 useFoPSF = True,      f0SD = 1,
                 useRGradient = False, maxGradFreq = 200,
                 noiseSD = 0.0,
                 showTissueSpectra = False):
        self.useFoPSF = useFoPSF                    # Use Off-Resonance Point Spread Function (PSF) to create Gaussian Fo Distributions
        self.useRGradient = useRGradient            # Use Linear Off-Resonance Gradient
        self.f0SD = f0SD                            # Off-Resonance Std. Dev. in Hz
        self.maxGradFreq = maxGradFreq              # Max Frequency for Off-Resonance Linear Gradient
        self.noiseSD = noiseSD                      # Signal Nosie Std. Dev.
        
        self.showTissueSpectra = showTissueSpectra  # Show Tissue Spectra
        pass

    def __str__(self):
        return "TissueOptions: \n UseFoPSF = %s,     FoSD = %f \n UseGradient = %s, MaxGradFreq = %f \n NoiseSD = %f \n ShowTissueSpectra = %s" % (self.useFoPSF, self.f0SD, self.useRGradient, self.maxGradFreq, self.noiseSD, self.showTissueSpectra)

class Tissue:
    def __init__(self, fileIndex = 0, tissueOptions = TissueOptions()):
        files = ['RectPhantom64Tissues.png', 'Phantom64Tissues.png', 'Phantom256Tissues.png']
        self.tissueTypes = {'Empty' : 0, 'Tissue1' : 255, 'Tissue2' : 91 }
        
        self.tissueOptions = tissueOptions
        self.loadTissueMask(files[fileIndex])
        self.generateTissue()
        self.generateResonancePSF(self.tissueOptions)
    
    def loadTissueMask(self, filename):
        imageMask = Image.open(filename).convert("L")
        imageMask = asmatrix(imageMask)
        self.Nx = shape(imageMask)[0]
        self.Ny = shape(imageMask)[1]
        
        self.mask = empty_like(imageMask)
        for r in range(self.Ny):
            for c in range(self.Nx):
                self.mask[r,c] = imageMask[r,c]
    
    def showTissueMask(self):
        #plt.imshow(self.mask)
        plt.imshow(self.mask, cmap = cm.Greys_r)
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
                self.tissue[r][c].f0 = self.tissue[r][c].f0 + tissueOptions.f0SD * random.uniform()
    
    def __str__(self):
        return "Tissue: \n Nx = %f, Ny = %f" % (self.Nx, self.Ny)

class PulseParameters:
    def __init__(self, downSampleFactor = 1.0,
                 FOVx = 0.256,     FOVy = 0.256,
                 dt = 1e-3,        tau = 1e-3,
                 B0 = 3.0,         gamma =  3 * 42.58e6,
                 TR = 10.0/1000.0, TE = 5.0/1000.0,
                 Nr = 100,
                 alpha = math.pi/5.0, phi = 0, dphi = math.pi/2.0):
        self.downSample = downSampleFactor        # Down Sampling Factor
        self.FOVx = FOVx                          # Field of View (FOV) in meters
        self.FOVy = FOVy                          # Field of View (FOV) in meters
        self.dt = dt                              # dt - Sample Period
        self.tau = tau                            # tau - Time for Phase Encode
        self.B0 = B0                              # B0 Field (T)
        self.gamma = gamma                        # Gyromagentic Ratio (Hz/T)
        self.TR = TR                              # TR - Time between Flips
        self.TE = TE                              # TE - Time from Flip to Sample
        self.Nr = Nr                              # Number of Repetitions
        self.alpha = alpha                        # Tip Angle
        self.phi = phi                            # Phi
        self.dphi = dphi                          # dPhi
        #self.generateParameters()                 # x, y, dx, dy, dkx, dky, kxMax, kyMax, GxMax, Gyi, GyMax
    
    '''
    def generateParameters(self):
        # Calculate Resolution and Position
        self.x = linspace(0, self.FOVx, self.Nx);
        self.y = linspace(0, self.FOVy, self.Ny);
        self.dx = self.FOVx / self.Nx;
        self.dy = self.FOVy / self.Ny;
    
        # Calculate kx and ky limits
        self.dkx = 1 / self.FOVx;
        self.dky = 1 / self.FOVy;
        self.kxMax = 1 / (2 * self.dx);
        self.kyMax = 1 / (2 * self.dy);
    
        # Calculate Gradient Amplitudes
        self.GxMax = 1/ ( (self.Gamma/(2*pi)) * self.FOVx * self.dt);
        self.Gyi = 1 / ( (self.Gamma/(2*pi)) * self.FOVy * self.tau);
        self.GyMax = self.Gyi * self.Ny / 2;
    '''
    
    def __str__(self):
        return "Pulse Parameters: \n DownSampleFactor=%f \n FOXx=%f,  FOXy=%f \n dt=%f,    tau=%f \n B0=%f,    gamma=%f \n Nr=%i,         TR=%f, TE=%f \n alpha=%f, phi=%f, dphi=%f" % (self.downSample, self.FOVx, self.FOVy, self.dt, self.tau, self.B0, self.gamma, self.Nr, self.TR, self.TE, self.alpha, self.phi, self.dphi)

class MRISimulator:

    def __init__(self):
        pass

    def generateMxySignal(self, tissue, sequenceParameters):

        # Get Sequence Parameters
        P = sequenceParameters;
        TR = P.TR;              # TR
        TE = P.TE;              # TE
        Nr = P.Nr;              # Number of Repetitions
        alpha = P.alpha;        # Tip Angle
        phi = P.phi;            # Phi
        dphi = P.dphi;          # dPhi

        
        Rx = mat([[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                  [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                  [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])
        
        # Get Tissue Parameters
        M0 = 1;
        T1 = tissue.T1;
        T2 = tissue.T2;
        b = 2 * pi * tissue.f0 * TR;

        # Generate Singal
        M = mat([0, 0, M0]).T
        for n in range(Nr): # Interate through Tips
            
            # Alpha Degree Tip
            Mtip = Rx * M
        
            # T1, T2 Relaxation
            E = mat(diag(array([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)])))
            Ez = mat([0,0,M0*(1-exp(-TR/T1))]).T
        
            # Off Resonance Precession
            P = mat([[cos(b),sin(b),0], [-sin(b),cos(b),0], [0,0,1]])
        
            # Singal after Relaxation and Off-Resonace Precession
            M = P * E * Mtip + Ez;
    
            # Increasing Phase for each exitation for SSPF 
            phi = phi + dphi
            Rx = mat([[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                      [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                      [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])
        
        # After one last tip, and T1, T2 Relaxation and Beta Precession
        Mtip = Rx * M;
        E = mat(diag(array([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)])))
        Ez = mat([0,0,M0*(1-exp(-TE/T1))]).T
        P = mat([[cos(b * TE/TR),sin(b * TE/TR),0], [-sin(b * TE/TR),cos(b * TE/TR),0], [0,0,1]])
        M = P * E * Mtip + Ez;
    
        # Save Sample Point after Steady State is Reached
        Mxy = array([M[0],M[1]])
        Mc =  complex(M[0], M[1])

        return Mc

    def generateMxyImage(self, tissue, sequenceParameters):
        
        # Generate Mxy signal for each sample location without Gradients
        Mc = zeros((tissue.Ny, tissue.Nx),'complex')
        count = 0
        for r in range(tissue.Nx):
            for c in range(tissue.Ny):
                count = count + 1
                if count % 100 == 0:
                    print count
                if tissue.tissue[r][c].Name == 'Empty':
                    Mc[r,c] = 0
                else:
                    Mc[r,c] = self.generateMxySignal(tissue.tissue[r][c], sequenceParameters)
        return Mc

    def twoPointDixon(self, image1, image2):
        water = 0.5 * (image1 + image2)
        fat = 0.5 * (image1 - image2)
        return water, fat
    
    def threePointDixon(self, image1, image2, image3):
        import numpy as np
        phi = angle(conj(image1) * image3) / 2.0;
        water = 0.5 * (image1 + image2 * np.exp(-1j * phi))
        fat = 0.5 * (image1 - image2 * np.exp(-1j * phi))
        return water, fat

    def plotMagnitudePhaseImage(self, image):
        mag = absolute(image).astype('float')
        phase = angle(image).astype('float')
        
        plt.subplot(211)
        plt.imshow(mag, cmap = cm.Greys_r)
        plt.axis('off')
        plt.subplot(212)
        plt.imshow(phase, cmap = cm.Greys_r)
        plt.axis('off')
        plt.show()

    def surfMagnitudePhaseImage(self, image):
        import numpy as np
        
        X = np.linspace(-1, 1, 64)
        Y = np.linspace(-1, 1, 64)
        X, Y = np.meshgrid(X, Y)
        mag = absolute(image).astype('float')
        phase = angle(image).astype('float')

        
        fig = plt.figure()
        ax = fig.add_subplot(211, projection='3d')
        Z = mag
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                               linewidth=0, antialiased=True)
        ax.set_zlim3d(np.min(Z), np.max(Z))
        fig.colorbar(surf)
        ax = fig.add_subplot(212, projection='3d')
        Z = phase
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                               linewidth=0, antialiased=True)
        ax.set_zlim3d(np.min(Z), np.max(Z))
        fig.colorbar(surf)
        plt.show()


if __name__ == '__main__':
    
    options = TissueOptions()
    tissue = Tissue(1, options)
    elements = TissueElements()
    df = abs(elements.getTissueElement('WhiteMatter').f0 -  elements.getTissueElement('Fat').f0)
    p = PulseParameters()
    p.TE = p.TE - 1 / (2 * df)
    p2 = PulseParameters()
    p3 = PulseParameters()
    p3.TE = p3.TE + 1 / (2 * df)

    '''
    print p
    print options
    print tissue
    print elements.getTissueElement(elements.tissueTypes[0])
    print elements.getTissueElement(elements.tissueTypes[1])
    print elements.getTissueElement(elements.tissueTypes[2])
    tissue.showTissueMask()
    '''

    sim = MRISimulator()
    image1 = sim.generateMxyImage(tissue, p)
    image2 = sim.generateMxyImage(tissue, p2)
    image3 = sim.generateMxyImage(tissue, p3)
    
    #sim.plotMagnitudePhaseImage(image)
    #sim.surfMagnitudePhaseImage(image)
    sim.surfMagnitudePhaseImage(image1)
    sim.surfMagnitudePhaseImage(image2)
    sim.surfMagnitudePhaseImage(image3)

    water, fat = sim.twoPointDixon(image1,image2)
    sim.plotMagnitudePhaseImage(water)
    sim.plotMagnitudePhaseImage(fat)
    
    water, fat = sim.threePointDixon(image1, image2, image3)
    sim.plotMagnitudePhaseImage(water)
    sim.plotMagnitudePhaseImage(fat)

