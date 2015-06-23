from __future__ import division
import cython
import numpy as np
cimport numpy as np
import math
import cmath

cdef extern from "math.h":
    double cos(double x)
    double sin(double x)
    double exp(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[np.double_t, ndim=1] SSPF_MxySignal(double M0, double alpha, double phi, double dphi, int Nr,
                                                    double TR, double TE, double T1, double T2, double f0):

    cdef float b = 2 * math.pi *f0 * TR
    cdef np.ndarray[np.double_t, ndim=2] Rx
    cdef np.ndarray[np.double_t, ndim=1] M = np.array([0,0,M0])
    cdef np.ndarray[np.double_t, ndim=1] Mtip = np.array([0.0,0.0,0.0])
    cdef np.ndarray[np.double_t, ndim=1] z = np.array([0.0, 0.0])

    Rx = np.array( [[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                    [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                    [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])

    for n in xrange(Nr):
        # Alpha Degree Tip
        #Mtip = Rx * M
        Mtip[0] = Rx[0,0] * M[0] + Rx[0,1] * M[1] + Rx[0,2] * M[2]
        Mtip[1] = Rx[1,0] * M[0] + Rx[1,1] * M[1] + Rx[1,2] * M[2]
        Mtip[2] = Rx[2,0] * M[0] + Rx[2,1] * M[1] + Rx[2,2] * M[2]


        # T1, T2 Relaxation
        M[0] = Mtip[0] * exp(-TR/T2)
        M[1] = Mtip[1] * exp(-TR/T2)
        M[2] = 1 + (Mtip[2] - 1) * exp(-TR/T1)

        # Off Resonace Precession
        # Complex Math - exp x = cos x + i sin x = [cos[x], sin[x]]
        # z = M * exp(1j * b) = [M[0], M[1]] * [cos b, sin b]
        z[0] = M[0] * cos(b) - M[1] * sin(b)
        z[1] = M[0] * sin(b) + M[1] * cos(b)
        M[0] = z[0]
        M[1] = z[1]

        # Increasing Phase for each excitation for SSPF
        phi = phi + dphi
        Rx[0,0] = cos(alpha)*sin(phi)**2+cos(phi)**2
        Rx[0,1] = (1-cos(alpha))*cos(phi)*sin(phi)
        Rx[0,2] = -sin(alpha)*sin(phi)
        Rx[1,0] = (1-cos(alpha))*cos(phi)*sin(phi)
        Rx[1,1] = cos(alpha)*cos(phi)**2+sin(phi)**2
        Rx[1,2] = sin(alpha)*cos(phi)
        Rx[2,0] = sin(alpha)*sin(phi)
        Rx[2,1] = -sin(alpha)*cos(phi)
        Rx[2,2] = cos(alpha)

    # After one last tip, and T1, T2 Relaxation and Beta Precession
    #Mtip = Rx * M
    Mtip[0] = Rx[0,0] * M[0] + Rx[0,1] * M[1] + Rx[0,2] * M[2]
    Mtip[1] = Rx[1,0] * M[0] + Rx[1,1] * M[1] + Rx[1,2] * M[2]
    Mtip[2] = Rx[2,0] * M[0] + Rx[2,1] * M[1] + Rx[2,2] * M[2]

    # return Mtip * exp(1j * b * TE / TR) = [Mtip[0], Mtip[1]] * exp(-TE/T2) * [cos b * TE / TR, sin b * TE / TR]
    M[0] = Mtip[0] * exp(-TE/T2) * cos(b * TE / TR) - Mtip[1] * exp(-TE/T2) * sin(b * TE / TR)
    M[1] = Mtip[0] * exp(-TE/T2) * sin(b * TE / TR) + Mtip[1] * exp(-TE/T2) * cos(b * TE / TR)
    return M

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[np.double_t, ndim=1] SSPF_MxySignal2(double M0, double alpha, double phi, double dphi, int Nr,
                                                    double TR, double TE, double T1, double T2, double f0):

    cdef float b = 2 * math.pi *f0 * TR
    cdef np.ndarray[np.double_t, ndim=2] Rx
    cdef np.ndarray[np.double_t, ndim=1] M = np.array([0, 0, M0])
    cdef np.ndarray[np.double_t, ndim=1] Mtip = np.array([0.0, 0.0, 0.0])
    cdef np.ndarray[np.double_t, ndim=2] E
    cdef np.ndarray[np.double_t, ndim=2] P

    Rx = np.array( [[cos(alpha)*sin(phi)**2+cos(phi)**2, (1-cos(alpha))*cos(phi)*sin(phi),   -sin(alpha)*sin(phi)],
                    [(1-cos(alpha))*cos(phi)*sin(phi),   cos(alpha)*cos(phi)**2+sin(phi)**2, sin(alpha)*cos(phi)],
                    [sin(alpha)*sin(phi),                -sin(alpha)*cos(phi),               cos(alpha)]])
    E = np.array([[exp(-TR/T2), 0, 0], [0, exp(-TR/T2), 0], [0, 0, exp(-TR/T1)]])
    P = np.array([[cos(b), sin(b), 0], [-sin(b), cos(b), 0], [0, 0, 1]])


    for n in xrange(Nr):
        # Alpha Degree Tip
        # Mtip = Rx.dot(M)
        Mtip[0] = Rx[0,0] * M[0] + Rx[0,1] * M[1] + Rx[0,2] * M[2]
        Mtip[1] = Rx[1,0] * M[0] + Rx[1,1] * M[1] + Rx[1,2] * M[2]
        Mtip[2] = Rx[2,0] * M[0] + Rx[2,1] * M[1] + Rx[2,2] * M[2]

        # Signal after Relaxation and Off-Resonace Precession
        # M = P * E * Mtip + Ez;
        M[0] = P[0,0] * E[0,0] * Mtip[0] + P[0,1] * E[1,1] * Mtip[1] + P[0,2] * E[2,2] * Mtip[2]
        M[1] = P[1,0] * E[0,0] * Mtip[0] + P[1,1] * E[1,1] * Mtip[1] + P[1,2] * E[2,2] * Mtip[2]
        M[2] = P[2,0] * E[0,0] * Mtip[0] + P[2,1] * E[1,1] * Mtip[1] + P[2,2] * E[2,2] * Mtip[2] + M0 * (1 - exp(-TR / T1))

        # Increasing Phase for each exitation for SSPF
        phi = phi + dphi
        Rx[0,0] = cos(alpha)*sin(phi)**2+cos(phi)**2
        Rx[0,1] = (1-cos(alpha))*cos(phi)*sin(phi)
        Rx[0,2] = -sin(alpha)*sin(phi)
        Rx[1,0] = (1-cos(alpha))*cos(phi)*sin(phi)
        Rx[1,1] = cos(alpha)*cos(phi)**2+sin(phi)**2
        Rx[1,2] = sin(alpha)*cos(phi)
        Rx[2,0] = sin(alpha)*sin(phi)
        Rx[2,1] = -sin(alpha)*cos(phi)
        Rx[2,2] = cos(alpha)

    # After one last tip, and T1, T2 Relaxation and Beta Precession
    # Mtip = Rx.dot(M);
    Mtip[0] = Rx[0,0] * M[0] + Rx[0,1] * M[1] + Rx[0,2] * M[2]
    Mtip[1] = Rx[1,0] * M[0] + Rx[1,1] * M[1] + Rx[1,2] * M[2]
    Mtip[2] = Rx[2,0] * M[0] + Rx[2,1] * M[1] + Rx[2,2] * M[2]

    E = np.array([[exp(-TE/T2), 0, 0], [0, exp(-TE/T2), 0], [0, 0, exp(-TE/T1)]])
    b = b * TE / TR
    P = np.array([[cos(b), sin(b), 0], [-sin(b), cos(b), 0], [0, 0, 1]])

    # M = P * E * Mtip + Ez;
    M[0] = P[0,0] * E[0,0] * Mtip[0] + P[0,1] * E[1,1] * Mtip[1] + P[0,2] * E[2,2] * Mtip[2]
    M[1] = P[1,0] * E[0,0] * Mtip[0] + P[1,1] * E[1,1] * Mtip[1] + P[1,2] * E[2,2] * Mtip[2]
    M[2] = P[2,0] * E[0,0] * Mtip[0] + P[2,1] * E[1,1] * Mtip[1] + P[2,2] * E[2,2] * Mtip[2] + M0 * (1 - exp(-TE / T1))

    return M

@cython.boundscheck(False)
@cython.wraparound(False)
def SSPF_OffResonanceSpectrum(double M0, double alpha, double phi, double dphi, int Nr,
                              double TR, double TE, double T1, double T2, int Ns, double BetaMax, double f0):

    print "---------------------------------------------------------------------------"
    print "Running SSPF Sequence with alpha = %f, dphi = %f, T1 = %f, T2 = %f, TR = %f" % (alpha, dphi, T1, T2, TR)

    cdef np.ndarray[np.double_t, ndim=1] beta = np.linspace(-BetaMax, BetaMax, Ns)
    cdef np.ndarray[np.double_t, ndim=1] f = beta / TR / (2 * math.pi)
    cdef np.ndarray[np.double_t, ndim=1] M
    cdef np.ndarray Mc = np.zeros(Ns,'complex')

    for n in range(Ns): # Interate through Beta values
           M = SSPF_MxySignal(M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, f[n] + f0)
           Mc[n] = complex(M[0], M[1])

    return f, Mc

@cython.boundscheck(False)
@cython.wraparound(False)
def SSFP_MxyImage(tissue, sequenceParameters):

    S = sequenceParameters
    cdef double alpha = S.alpha
    cdef double phi = S.phi
    cdef double dphi = S.dphi
    cdef int Nr = S.Nr
    cdef double TR = S.TR
    cdef double TE = S.TE
    cdef double T1, T2, M0, f0
    cdef np.ndarray[np.double_t, ndim=1] M


     # Generate Mxy signal for each sample location without Gradients
    cdef np.ndarray Mc = np.zeros((tissue.Ny, tissue.Nx),'complex')
    for r in range(tissue.Nx):
        for c in range(tissue.Ny):
            if tissue.tissue[r][c].Name == 'Empty':
                Mc[r, c] = 0
            else:
                 T1 = tissue.tissue[r][c].T1
                 T2 = tissue.tissue[r][c].T2
                 M0 = tissue.tissue[r][c].protonDensity
                 f0 = tissue.tissue[r][c].f0
                 M =  SSPF_MxySignal(M0, alpha, phi, dphi, Nr, TR, TE, T1, T2, f0)
                 Mc[r, c] = complex(M[0], M[1])

    return Mc

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[np.double_t, ndim=1] ComplexMult(double real1, double imag1, double real2, double imag2):
    cdef np.ndarray[np.double_t, ndim=1] c = np.array([0.0, 0.0])
    c[0] = real1 * real2 - imag1 * imag2
    c[1] = real1 * imag2 + real2 * imag1
    return c

def ComplexMultTest(double real1, double imag1, double real2, double imag2):
    print real1, imag1, real2, imag2
    cdef np.ndarray[np.double_t, ndim=1] c = ComplexMult(real1,imag1,real2,imag2)
    print[c[0],c[1]]