import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("prueba.txt")

plt.figure()
plt.loglog(data[:,0], data[:,3], 'o-')
plt.show()