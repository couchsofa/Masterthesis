import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.legend_handler import HandlerLine2D

def AWS(ag, gamma1, n, S, TB, TC, TD, T):

	if T == 0:
		Se = ag * gamma1 * S * n * 2.5

	elif T <= TB:
		Se = ag * gamma1 * S * (1 + T/TB * (n * 2.5 -1))

	elif T <= TC:
		Se = ag * gamma1 * S * 2.5 * n

	elif T <= TD:
		Se = ag * gamma1 * S * n * 2.5 * TC/T

	else:
		Se = ag * gamma1 * S * n * 2.5 * (TC*TD)/(T*T)

	return Se

ag = 3.924
gamma1 = 1

S  = 1
TB = 0.4
TC = 1.6
TD = 2.0

_T = np.arange(0.01, 5.1, 0.01)
Sd = []
Sd_xi = []

for T in _T:
	n = 1
	Sd_ = AWS(ag, gamma1, n, S, TB, TC, TD, T)
	Sd.append(Sd_)

	n = 0.723
	Sd_ = AWS(ag, gamma1, n, S, TB, TC, TD, T)
	Sd_xi.append(Sd_)

	

plt.plot(_T, Sd, label=r'$S_a (\xi = 5\%)$', linestyle='-')
plt.plot(_T, Sd_xi, label=r'$S_a (\xi = 14.147\%)$', linestyle='--')

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=2, frameon=False)

plt.xlabel('T [$s$]')
plt.ylabel('$S_a$ $[m/s^{2}]$')

plt.grid(True)
plt.tight_layout()
#plt.axes().set_aspect(0.1)
plt.gca().set_ylim(0,11)
plt.margins(x=0, y=0)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/AWS_beispiel.png", dpi=300)
plt.clf()
