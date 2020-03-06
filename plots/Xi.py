import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.legend_handler import HandlerLine2D

_Xi = np.arange(0.001, 25.01, 0.1)
#_T = [1]

Norm = []
for Xi in _Xi:
	Norm.append(np.sqrt(10/(5+Xi)))
	
V = []
for Xi in _Xi:
	V.append(10/(2*Xi))

VT = []
for Xi in _Xi:
	Xi = Xi /100
	VT.append(0.1/(2*Xi*np.sqrt(1-Xi)))	

plt.plot(_Xi, Norm, label=r'$\eta = \sqrt{\frac{10}{5+\xi}}$', linestyle=':')
plt.plot(_Xi, V, label=r'$\eta = \frac{10}{2\xi}$', linestyle='-')
plt.plot(_Xi, VT, label=r'$\eta = \frac{1}{2\frac{\xi}{10}\sqrt{1-\frac{\xi}{100}}}$', linestyle='--')


plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=3, frameon=False)

plt.xlabel(r'$\xi  [\%]$')
plt.ylabel(r'$\eta  [-]$')

plt.grid(True)
plt.tight_layout()
#plt.axes().set_aspect(0.7)
plt.margins(x=0, y=0.1)
plt.gca().set_ylim(0,3)
plt.gca().set_xlim(0,25)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/Xi.png", dpi=300)
plt.clf()