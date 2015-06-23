import numpy as np
import matplotlib.cm as cm
from matplotlib.pylab import plt
from math import *

def Test():
    list = [1,2,3]
    a = np.array(list)
    b = np.array([[1,2,3], [4,5,6], [7,8,9]])
    c1 = b * b       # .* - MATLAB
    c2 = b.dot(a)    # * - MATLAB
    print c1
    print c2

def Test2():
    t = np.linspace(0,10,100)
    x = np.sin(t)

    plt.plot(t, x)
    plt.title('sin(x)')
    plt.xlabel('time (s)')
    plt.ylabel('magnitude')
    plt.show()

def MRI():

    T1 = 800*10**-3
    T2 = 100*10**-3
    a = pi/2
    t = np.linspace(0, 1, 100)

    # Initalize and Rotate
    M = np.array([0,0,1])
    R = np.array([[1, 0,       0],
                  [0, cos(a),  sin(a)],
                  [0, -sin(a), cos(a)]])
    M = R.dot(M)

    # T1/T2 Relaxationq
    M_z = 1 + (M[2] - 1) * np.exp(-t/T1)
    M_x = M[0]*np.exp(-t/T2)
    M_y = M[1]*np.exp(-t/T2)
    M_xy = np.sqrt(M_x*M_x + M_y*M_y)

    plt.subplot(121)
    plt.plot(t, M_xy)
    plt.subplot(122)
    plt.plot(t, M_z)
    plt.show()

def SSFP():
    T1 = 800*10**-3
    T2 = 100*10**-3
    a = pi/2
    Tr = .01
    t = np.linspace(0, Tr, 100)
    M = np.array([0,0,1])
    R = np.array([[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]])

    Mx = []
    My = []
    Mz = []

    Nr = 10
    for n in range(Nr):

        # Initalize and Rotate
        M = R.dot(M)
        # T1/T2 Relaxationq
        M_z = 1 + (M[2] - 1) * np.exp(-t/T1)
        M_x = M[0]*np.exp(-t/T2)
        M_y = M[1]*np.exp(-t/T2)
        M_xy = np.sqrt(M_x*M_x + M_y*M_y)
        M[0] = M_x[-1]
        M[1] = M_y[-1]
        M[2] = M_z[-1]

        Mx = np.append(Mx, M_x)
        My = np.append(My, M_y)
        Mz = np.append(Mz, M_z)

    tnew = np.linspace(0, 10*Tr, 1000)
    plt.subplot(131)
    plt.plot(tnew, Mx)
    plt.subplot(132)
    plt.plot(tnew, My)
    plt.subplot(133)
    plt.plot(tnew, Mz)
    plt.show()


if __name__ == '__main__':
    SSFP()

