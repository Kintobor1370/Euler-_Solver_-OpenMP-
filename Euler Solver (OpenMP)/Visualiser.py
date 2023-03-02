import matplotlib.pyplot as plt
import numpy as np
from array import array

file = open('Data\VData.txt', 'r')
ExperimentNumber = float(file.readline())
xLabel = file.readline()
yLabel = file.readline()
size = int(file.readline())

x = array("d")
y = array("d")
eq = array("d")
x = np.fromfile(file, "d", size, "\n")
y = np.fromfile(file, "d", size, "\n")
eq = np.full(size, 200)

plt.plot(x, y)
plt.ylim(0, 400)
plt.xlabel(xLabel)
plt.ylabel(yLabel)
if ExperimentNumber == 14.6:
    plt.plot(x, eq, linestyle = '--', color = 'grey')

print("\n============================PLOT CREATED!=============================\n")
plt.show()
