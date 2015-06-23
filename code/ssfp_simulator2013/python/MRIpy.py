# ------------ Tissue (ms) -----------
# Gray Matter  : T1 - 920, T2 - 100
# White Matter : T1 - 790, T2 - 92
# Fat:           T1 - 270, T2 - 85
# Muscle:        T1 - 870, T2 - 47
# ------------ Tissue (ms) ------------
# Fat Water Spin Separation = 228 Hz

from math import *
from cmath import *
from numpy import *

import matplotlib.cm as cm
from matplotlib.pylab import plt

import MRI
import Tissue



class MRISimulator():
    
    def __init__(self):
        pass

    useOptimized = True
    modes = 'iterative', 'steadystate', 'iterativeMatrix'
    mode = 'iterative'
    debug = True
    
    def SSFP_MxySignal(self, M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, f0):
        
        Rx = mat([[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                  [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                  [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])
        
        # Get Off-Resonace Phase 
        b = 2 * pi * f0 * TR
        
        if MRISimulator.mode == 'iterative':
            M = mat([0, 0, M0]).T
            for n in xrange(Nr):
                
                # Alpha Degree Tip
                Mtip = Rx * M
                
                # T1, T2 Relaxation
                M = matrix([ Mtip[0,0] * exp(-TR/T2),
                            Mtip[1,0] * exp(-TR/T2),
                            1 + (Mtip[2,0] - 1) * exp(-TR/T1)]).T

                # Off Resonace Precession
                z = complex(M[0], M[1])
                z = z * exp(1j * b)
                M[0] = real(z)
                M[1] = imag(z)
            
                # Increasing Phase for each excitation for SSPF
                phi = phi + dphi
                Rx = mat([[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                          [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                          [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])

            # After one last tip, and T1, T2 Relaxation and Beta Precession
            Mtip = Rx * M
            z = complex(Mtip[0,0] * exp(-TE/T2),
                        Mtip[1,0] * exp(-TE/T2))
            Mc = z * exp(1j * b * TE/TR)

        elif MRISimulator.mode == 'iterativeMatrix':
            # Generate Singal
            M = mat([0, 0, M0]).T
            for n in xrange(Nr): # Interate through Tips
            
                # Alpha Degree Tip
                Mtip = Rx * M
            
                # T1, T2 Relaxation
                E = mat(diag(array([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)])))
                Ez = mat([0,0,M0*(1-exp(-TR/T1))]).T
            
                # Off Resonance Precession
                P = mat([[cos(b),sin(b),0], [-sin(b),cos(b),0], [0,0,1]])
            
                # Signal after Relaxation and Off-Resonace Precession
                M = P * E * Mtip + Ez
            
                # Increasing Phase for each exitation for SSPF
                phi = phi + dphi
                Rx = mat([[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                          [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                          [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])
        
            # After one last tip, and T1, T2 Relaxation and Beta Precession
            Mtip = Rx * M
            E = mat(diag(array([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)])))
            Ez = mat([0,0,M0*(1-exp(-TE/T1))]).T
            P = mat([[cos(b * TE/TR),sin(b * TE/TR),0], [-sin(b * TE/TR),cos(b * TE/TR),0], [0,0,1]])
            M = P * E * Mtip + Ez
    
            # Save Sample Point after Steady State is Reached
            Mxy = array([M[0],M[1]])
            Mc =  complex(M[0], M[1])
    
        elif MRISimulator.mode == 'steadystate':
            # Generate Singal
            M = mat([0, 0, M0]).T
            I = mat([[1,0,0],[0,1,0],[0,0,1]])
            # T1, T2 Relaxation
            E = mat(diag(array([exp(-TR/T2), exp(-TR/T2), exp(-TR/T1)])))
            Ete = mat(diag(array([exp(-TE/T2), exp(-TE/T2), exp(-TE/T1)])))
            # Off Resonance Precesion
            P = mat([[cos(b),sin(b),0], [-sin(b),cos(b),0], [0,0,1]])
            Pte = mat([[cos(b*TE/TR),sin(b*TE/TR),0], [-sin(b*TE/TR),cos(b*TE/TR),0], [0,0,1]])
    
            Mminus = linalg.inv(I - P * E * Rx) * (I - E) * M        # Before Last Tip
            Mplus = Rx * Mminus                                      # Singal After Tip
            Mte = Pte * Ete * Mplus + (I - Ete) * M                  # After Last Relaxation and Precession
            M = Mte
    
            # Save Sample Point after Steady State is Reached
            Mxy = array([M[0],M[1]])
            Mc =  complex(M[0], M[1])
                
        return Mc
    
    def SSFP_OffResonanceSpectrum(self, M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, Ns, BetaMax, f0):
    
        if MRISimulator.debug:
            print "---------------------------------------------------------------------------"
            print "Running SSPF Sequence with alpha = %f, dphi = %f, T1 = %f, T2 = %f, TR = %f" % (alpha, dphi, T1, T2, TR)

        if MRISimulator.useOptimized:
            return MRI.SSPF_OffResonanceSpectrum(M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, Ns, BetaMax, f0)
        else:
            beta = linspace(-BetaMax, BetaMax, Ns)
            f = beta / TR / (2 * pi)
            Mc = zeros(Ns, 'complex')

            for n in range(Ns): # Interate through Beta values
                Mc[n] = self.SSFP_MxySignal(M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, f[n] + f0)

            return f, Mc

    def SSFP_MxyImage(self, tissue, sequenceParameters):


        if MRISimulator.useOptimized:
            return MRI.SSFP_MxyImage(tissue, sequenceParameters)

        S = sequenceParameters
        alpha = S.alpha
        phi = S.phi
        dphi = S.dphi
        Nr = S.Nr
        TR = S.TR
        TE = S.TE

        # Generate Mxy signal for each sample location without Gradients
        Mc = zeros((tissue.Ny, tissue.Nx),'complex')
        for r in range(tissue.Nx):
            for c in range(tissue.Ny):
                if tissue.tissue[r][c].Name == 'Empty':
                    Mc[r, c] = 0
                else:
                    T1 = tissue.tissue[r][c].T1
                    T2 = tissue.tissue[r][c].T2
                    M0 = tissue.tissue[r][c].protonDensity
                    f0 = tissue.tissue[r][c].f0
                    Mc[r, c] = self.SSFP_MxySignal(M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, f0)
        return Mc

    def showMagPhasePlot(self, t, x):
        mag = absolute(x)
        phase = angle(x)
    
        plt.subplot(211)
        plt.plot(t, mag)
        plt.ylabel('Magitude')
        plt.title('SSPF Sequence')
        plt.grid(True)

        plt.subplot(212)
        plt.plot(t, phase)
        plt.xlabel('Off-Resonance (Hz)')
        plt.ylabel('Phase')
        plt.grid(True)
        plt.show()

    def saveMagPhasePlot(self, t, x, filename):
        mag = absolute(x).astype('float')
        phase = angle(x).astype('float')
        
        plt.subplot(211)
        plt.plot(t, mag)
        plt.ylabel('Magitude')
        plt.title('SSPF Sequence')
        plt.grid(True)
        
        plt.subplot(212)
        plt.plot(t, phase)
        plt.xlabel('Off-Resonance (Hz)')
        plt.ylabel('Phase')
        plt.grid(True)
        plt.savefig(filename)

    def showMagPhaseImage(self, image):
        mag = absolute(image).astype('float')
        phase = angle(image).astype('float')

        plt.subplot(211)
        plt.imshow(mag, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(212)
        plt.imshow(phase, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.show()

    def saveMagPhaseImage(self, image, filename):
        mag = absolute(image).astype('float')
        phase = angle(image).astype('float')

        plt.subplot(211)
        plt.imshow(mag, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(212)
        plt.imshow(phase, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.savefig(filename)

    def saveMagPhaseImage2(self, image1, image2, filename):
        mag = absolute(image1).astype('float')
        phase = angle(image1).astype('float')
        mag2 = absolute(image2).astype('float')
        phase2 = angle(image2).astype('float')

        plt.subplot(221)
        plt.imshow(mag, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(222)
        plt.imshow(phase, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.subplot(223)
        plt.imshow(mag2, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(224)
        plt.imshow(phase2, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.savefig(filename)

    def saveMagPhaseImage3(self, image1, image2, image3, filename):
        mag = absolute(image1).astype('float')
        phase = angle(image1).astype('float')
        mag2 = absolute(image2).astype('float')
        phase2 = angle(image2).astype('float')
        mag3 = absolute(image3).astype('float')
        phase3 = angle(image3).astype('float')

        plt.subplot(321)
        plt.imshow(mag, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(322)
        plt.imshow(phase, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.subplot(323)
        plt.imshow(mag2, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(324)
        plt.imshow(phase2, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.subplot(325)
        plt.imshow(mag3, cmap=cm.Greys_r)
        plt.title('Magnitude')
        plt.axis('off')
        plt.subplot(326)
        plt.imshow(phase3, cmap=cm.Greys_r)
        plt.title('Phase')
        plt.axis('off')
        plt.savefig(filename)

def SpectrumTest(useOptimized = False):
    T1 = 790.0/1000.0; T2 = 92.0/1000.0;
    T1F = 270.0/1000.0; T2F = 85.0/1000.0
    T1M = 870.0/1000.0; T2M = 47.0/1000.0

    sim = MRISimulator()
    MRISimulator.useOptimized = useOptimized
    f0, Mc = sim.SSFP_OffResonanceSpectrum(M0 = 1,
                                             alpha = math.pi/3.0,
                                             phi = 0.0,
                                             dphi= 0.0,
                                             Nr = 200,
                                             TR = 10.0/1000.0,
                                             TE = 5.0/1000.0,
                                             T1 = T1,
                                             T2 = T2,
                                             Ns = 200,
                                             BetaMax = math.pi,
                                             f0=0.0)
    sim.showMagPhasePlot(f0, Mc)

def ImageTest(useOptimized):
    options = Tissue.TissueOptions()
    tissue = Tissue.Tissue(1, options)
    elements = Tissue.TissueElements()
    p = Tissue.PulseParameters()
    sim = MRISimulator()
    MRISimulator.useOptimized = useOptimized
    img = sim.SSFP_MxyImage(tissue, p)
    sim.showMagPhaseImage(img)

if __name__ == '__main__':
    pass
