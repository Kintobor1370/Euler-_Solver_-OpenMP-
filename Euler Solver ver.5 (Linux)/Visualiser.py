import matplotlib.pyplot as plt
import matplotlib.patches as mp
import numpy as np
from array import array

file = open('Data/VisData.txt', 'r')
ExperimentNumber = float(file.readline())
title = file.readline()
xLabel = file.readline()
yLabel = file.readline()
size = int(file.readline())

x = array("d")
y1 = array("d")
y2 = array("d")
eq = array("d")
x = np.fromfile(file, "d", size, "\n")
y1 = np.fromfile(file, "d", size, "\n")
y2 = np.fromfile(file, "d", size, "\n")
eq = np.full(size, 200)

plt.plot(x, y1, color='blue')
plt.plot(x, y2, color='orange')
plt.title(title)
plt.xlabel(xLabel)
plt.ylabel(yLabel)
plt.ylim(0, 350)

if ExperimentNumber == 14.6:
    plt.xlim(0, 0.2)
    plt.ylim(0, 400)
    plt.plot(x, eq, linestyle = '--', color = 'grey')

plt.show()



er = [0] * size
for i in range(size):
    er[i] = abs(y1[i] - y2[i])
    er[i] = er[i] / y2[i]
    er[i] = er[i] * 100


plt.plot(x, er, color='blue')
plt.title(title)
plt.xlabel(xLabel)
plt.ylabel("Error (in %)")
plt.ylim(0, 10)

print("\n============================PLOTS CREATED!=============================\n")
plt.show()
