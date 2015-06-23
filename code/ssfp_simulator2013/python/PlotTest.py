#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import Image


if False:
    # Generate image
    N=25
    data=np.random.random_sample((N,N))
    plt.imshow(data, interpolation='nearest')
    plt.show()

elif False:
    # Load Image
    fname = 'Phantom256Tissues.png'
    image = Image.open(fname).convert("L")
    arr = np.asarray(image)
    
    plt.subplot(211)
    plt.imshow(arr, cmap = cm.Greys_r)
    plt.axis('off')
    #plt.colorbar()
    
    fname = 'Phantom64Tissues.png'
    image = Image.open(fname).convert("L")
    arr = np.asarray(image)
    
    plt.subplot(212)
    plt.imshow(arr, cmap = cm.Greys_r)
    plt.axis('off')
    #plt.colorbar()
    plt.show()

elif False:

    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.arange(-10, 10, 0.25)
    Y = np.arange(-10, 10, 0.25)
    X, Y = np.meshgrid(X, Y)
    R = np.sqrt(X**2 + Y**2)
    Z = np.sin(R)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=False)
    plt.show()

else:
    
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(-1, 1, 64)
    Y = np.linspace(-1, 1, 64)
    X, Y = np.meshgrid(X, Y)
    fname = 'Phantom64Tissues.png'
    image = Image.open(fname).convert("L")
    Z = np.asarray(image)

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=True)
    ax.set_zlim3d(np.min(Z), np.max(Z))
    fig.colorbar(surf)
    plt.show()


