#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib as mp
import sys
import math
import numpy as np
import locale

from mpl_toolkits.axes_grid1 import make_axes_locatable

def read_file(filename):
    f = open(filename, 'r')
    sizes = f.readline().split()
    nx = int(sizes[0])
    ny = int(sizes[1])
   
    ey = np.zeros((nx,ny))
    bx = np.zeros((nx,ny))
    ne = np.zeros((nx,ny))
    
    for i in range(nx*ny):
        values = f.readline().split()
        x = int(values[0])
        y = int(values[1])
        ey[x,y] = float(values[2])
        bx[x,y] = float(values[3])
        ne[x,y] = float(values[4])
        
    f.close()
    return ey.transpose(), bx.transpose(), ne.transpose()

def main():
    locale.setlocale(locale.LC_ALL, '')
    nmin = int(sys.argv[1])
    nmax = int(sys.argv[2])
    delta = int(sys.argv[3])

    for i in range(nmin, nmax+1, delta):
        ey, bx, ne = read_file('iteration_%d half.txt' % i)
        fig = plt.figure(num=None, figsize=(10,5))
              
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
        
        ax1.set_title("Ey")
        ax2.set_title("Ne")
                
        surf = ax1.imshow(-ey, interpolation = None, origin = 'lower', cmap = 'bwr')
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(surf,  orientation  = 'vertical', cax=cax)
        
        surf = ax2.imshow(ne, interpolation = None, origin = 'lower', cmap = 'jet')
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(surf,  orientation  = 'vertical', cax=cax)
        
        plt.tight_layout()
        plt.savefig('iteration_%d half.png' % i)
        plt.close()

if __name__ == '__main__':
    main()
