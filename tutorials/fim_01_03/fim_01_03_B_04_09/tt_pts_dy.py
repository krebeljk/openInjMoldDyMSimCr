import numpy as np
import sys

# stevilo tock
pts = 64

# debelina stene
L = 0.00121 + float(sys.argv[1])

z1 = L/pts/2
z2 = L - z1

myZ = np.linspace(z1, z2, pts)

def myLine(i):
     return "(37.9e-3 0.05e-3 " + str(myZ[i]) + ")"

myLines = []

for i in range(pts):
    myLines.append(myLine(i))

ttPts = "\n".join(myLines)

with open("tt_pts_dy.txt", "w") as txt:
    txt.write(ttPts)
