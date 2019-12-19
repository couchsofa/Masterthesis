import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.legend_handler import HandlerLine2D

def G(m, c, k, T):
	omega = (2 * math.pi) / T
	s = omega * 1j
	G = 1 / (m*s*s + c*s + k)
	return np.absolute(G)


def AWS(ag, gamma1, beta0, q, S, TB, TC, TD, T):

	if T == 0:
		Se = ag * gamma1 * S

	elif T <= TB:
		Se = ag * gamma1 * S * (1 + T/TB * (beta0/q -1))

	elif T <= TC:
		Se = ag * gamma1 * S * beta0/q

	elif T <= TD:
		Se = ag * gamma1 * S * beta0/q * TC/T

	else:
		Se = ag * gamma1 * S * beta0/q * (TC*TD)/(T*T)

	return Se


m = 0.2
k = 5.0
c = 0.05*2*math.sqrt(m*k)

ag = 0.4
gamma1 = 1
beta0 = 2.5
q = 1.5

S = 0.75
TB = 0.1
TC = 0.5
TD = 2

_T = np.arange(0.01, 3.0, 0.01)
Sd = []
G_abs = []
SdG = []

for T in _T:
	Sd_ = AWS(ag, gamma1, beta0, q, S, TB, TC, TD, T)
	G_ = G(m, c, k, T)

	Sd.append(Sd_)
	G_abs.append(G_)
	SdG.append(Sd_ * G_)

plt.plot(_T, Sd, label='$S_d$ $[m/s^{2}]$', linestyle='-')
plt.plot(_T, G_abs, label='$|G|$ $[-]$', linestyle=':')
plt.plot(_T, SdG, label='$S_d*|G|$ $[m/s^{2}]$', linestyle='--')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=3, frameon=False)

plt.xlabel('T [$s$]')

plt.grid(True)
plt.tight_layout()
plt.axes().set_aspect(0.7)
plt.margins(x=0, y=0)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/AWS.png", dpi=300)
plt.clf()
