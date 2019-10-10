import numpy as np
from ROOT import *
import sys
a = np.loadtxt(sys.argv[1])
NY = a.shape[1]
NX = a.shape[1]
NZ = a.shape[0]/a.shape[1]

LX = 24.00
LY = 24.00
LZ = 9.00
dx = LX/NX
dy = LY/NY
dz = LZ/NZ
Xmin,Xmax = -LX/2, LX/2
Ymin,Ymax = -LY/2, LY/2
Zmin,Zmax = -LZ/2, LZ/2
print Xmin+ NX*dx
print Ymin+ NY*dy
print Zmin+ NZ*dz
print NX,NY,NZ, dx, dy, dz
A = a.reshape(NZ,NY,NX)
A = A[1:-1]
NZ = NZ-2
A = np.rot90(A,2, axes= (1,2))
A[A<0.3] = 0
#for nb,i in enumerate(A):
#    plt.imshow(i)
#    print nb*dz + Zmin#
#    plt.show()

f = TFile('reconstruction.root','recreate')
RSPHist = TH3D('RSP','RSP',NX, Xmin, Xmax, NY, Ymin, Ymax, NZ, Zmin, Zmax)
b = np.zeros((NZ+2,NY+2,NX+2),'float64')
b[1:-1,1:-1,1:-1] = A

b = b.flatten()
RSPHist.Set(len(b),b)
RSPHist.Write()
f.Close()
