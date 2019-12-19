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


m2 = 1
k2 = 50

m1 = 0.2
#k1 = 5.0

ag = 0.4
gamma1 = 1
beta0 = 2.5
q = 1.5

S = 0.75
TB = 0.1
TC = 0.5
TD = 2

def c(m1, k1):
	return 0.05*2*math.sqrt(m1*k1)

# 2MS
def _2MS_F2(m1, k1, m2, ag, gamma1, beta0, q, S, TB, TC, TD):
	omega_1 = math.sqrt(k1/(m1+m2))
	T_1 = (2*math.pi)/omega_1
	Sd_1 = AWS(ag, gamma1, beta0, q, S, TB, TC, TD, T_1)
	Fb_1 = (m1+m2)*Sd_1
	X = (k1+k2-m1*omega_1*omega_1)/k2
	X_norm = X/math.sqrt(1+X*X)

	return [X_norm*Fb_1, omega_1]

# 1MS
def _1MS_F2(m1, k1, c, m2, ag, gamma1, beta0, q, S, TB, TC, TD):
	omega_2 = math.sqrt(k2/m2)
	T_2 = (2*math.pi)/omega_2
	Sd_2 = AWS(ag, gamma1, beta0, q, S, TB, TC, TD, T_2) * G(m1, c, k1, T_2)

	return [Sd_2, omega_2]



k1_k2 = np.arange(0.06, 0.15, 0.001)
_2MS_to_1MS = []
omega_1_to_omega_2 = []

for a in k1_k2:
	c_ = c(m1, k2*a)

	_2MS = _2MS_F2(m1, k2*a, m2, ag, gamma1, beta0, q, S, TB, TC, TD)[0]
	_1MS = _1MS_F2(m1, k2*a, c_, m2, ag, gamma1, beta0, q, S, TB, TC, TD)[0]

	_2MS_to_1MS.append(_1MS/_2MS)

	omega_1_to_omega_2.append(_2MS_F2(m1, k2*a, m2, ag, gamma1, beta0, q, S, TB, TC, TD)[1]/_1MS_F2(m1, k2*a, c_, m2, ag, gamma1, beta0, q, S, TB, TC, TD)[1])

plt.plot(omega_1_to_omega_2, _2MS_to_1MS)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=3, frameon=False)

plt.xlabel(r'$\omega1/\omega2 [-]$')
plt.ylabel('2MS / 1MS [$-$]')

plt.grid(True)
plt.tight_layout()
plt.axes().set_aspect(0.2)
plt.margins(x=0, y=0)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/DampenerSweep.png", dpi=300)
plt.clf()
