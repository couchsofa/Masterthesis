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
k2 = 55.8

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

def _c(k1):
	return 0.05#*2*math.sqrt(m1*k1)

# 2MS
def _2MS_F2(k1):
	#omega_1 = math.sqrt(k1/(m1+m2))

	omega_1 = math.sqrt( (m1*k2+m2*(k1+k2) - math.sqrt( (m1*k2+m2*(k1+k2))**2 - 4*m1*m2*( (k1+k2)*k2 - k2**2 ) )) / (2*m1*m2) )

	T_1 = (2*math.pi)/omega_1
	Sd_1 = AWS(ag, gamma1, beta0, q, S, TB, TC, TD, T_1)
	Fb_1 = m2*Sd_1
	X = (k1+k2-m1**omega_1)/k2
	X_norm = 1/X

	return X_norm*Fb_1

# 1MS
def _1MS_F2(k1):
	c = _c(k1)
	omega_2 = math.sqrt(k2/m2)
	T_2 = (2*math.pi)/omega_2
	Sd_2 = m2*AWS(ag, gamma1, beta0, q, S, TB, TC, TD, T_2) * G(m1, c, k1, T_2)

	return Sd_2



k1_k2 = np.arange(0.06, 0.15, 0.001)
_1MS = []
_2MS = []
_1MS_to_2MS = []
omega_1_to_omega_2 = []

for a in k1_k2:
	k1 = k2*a

	omega1 = math.sqrt(k1/(m2+m1))
	omega2 = math.sqrt(k2/m2)

	_2MS.append(_2MS_F2(k1))
	_1MS.append(_1MS_F2(k1))
	_1MS_to_2MS.append(_1MS_F2(k1)/_2MS_F2(k1))
	omega_1_to_omega_2.append(omega1/omega2)

plt.plot(omega_1_to_omega_2, _1MS_to_2MS, linestyle='-')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=3, frameon=False)

plt.xlabel(r'$\omega1/\omega2 [-]$')
plt.ylabel('1MS / 2MS [$-$]')

plt.grid(True)
plt.tight_layout()
plt.axes().set_aspect(0.2)
plt.margins(x=0, y=0)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/DampenerSweep.png", dpi=300)
plt.clf()

# Im Bereich der Annahme, dass das untere System dominiert 20% Fehler, wird die Struktur ggü der Basis weicher/steifer
# so überlagern sich die Eigenfrequenzen (Schwebungen, Angewandete Baudynamik S. 100)
# Plot für Dämpfung = 5%

# Resoanzbereiche da 2MS hier die ungedämpfte Schwingung betrachtet