import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.legend_handler import HandlerLine2D

#
#	m1
#	|	
#	|
#	m2
#	|
#	|
#  ---
#

# Elastisches Antwortspektrum
def Se(n, ag, S, TB, TC, TD, T):

	if T <= TB:
		Se = ag * S * (1 + T/TB * (n*2.5 -1))

	elif T <= TC:
		Se = ag * S * 2.5 * n

	elif T <= TD:
		Se = ag * S * n * 2.5 * TC/T

	else:
		Se = ag * S * n * 2.5 * (TC*TD)/(T*T)

	return Se


def vereinfacht_rayleigh(Xi1, k1, m1, Xi2, k2, m2):

	# Omega1 ungedämpft
	A = (k2 + k1) * m1 + k1 * m2
	B = ((k2 + k1) * m1 + k1 * m2)**2 - 4 * m2 * m1 * k2 * k1
	C = 2 * m2 * m1

	omega1 = np.sqrt((A - np.sqrt(B))/(C))
	omega2 = np.sqrt((A + np.sqrt(B))/(C))

	# Umrechnung mit Rayleigh Dämpfung
	e1 = (k2 + k1 - m2 * omega1**2)/k1
	phi_11 = np.sqrt(1/(1+e1**2))
	phi_21 = e1 * phi_11

	e2 = (k2 + k1 - m2 * omega2**2)/k1
	phi_12 = np.sqrt(1/(1+e2**2))
	phi_22 = e2 * phi_12

	_m2 = phi_11**2 * m2 + phi_21**2 * m1
	_k2 = phi_11**2 * (k2 + k1) - 2 * phi_21 * phi_11 * k1 + phi_21**2 * k1

	_omega1 = np.sqrt(_k2/_m2)

	_m1 = phi_12**2 * m2 + phi_22**2 * m1
	_k1 = phi_12**2 * (k2 + k1) - 2 * phi_22 * phi_12 * k1 + phi_22**2 * k1

	_omega2 = np.sqrt(_k1/_m1)

	# Eigenfrequenz erste Eigenform (Isolator) gedämpft
	omega_1d = _omega1 * np.sqrt(1 - Xi2**2)
	T1 = 2*np.pi/omega_1d

	# Transmissionskoeffizient des Systems

	# Nach Pocanschi
	a = (2*_omega1*_omega2*(Xi2 * _omega2 - Xi1 * _omega1))/(_omega2**2 - _omega1**2)
	b = (2*(Xi1 * _omega2 - Xi2 * _omega1))/(_omega2**2 - _omega1**2)

	#  Nach Isemann
	#b = 2 * Xi2 / (omega1 * omega2)
	#a = omega1 * omega2 * b

	_c1 = a * _m1 + b * _k1
	_c2 = a * _m2 + b * _k2
	

	VT = VT_1_over_3(omega_1d, _m2, _k2, _c2, _m1, _k1, _c1) 

	n = 1#np.sqrt(10/(5+Xi2*100))

	return Se(n, ag, S, TB, TC, TD, T1) * VT #/ 10.05


def VT_1_over_3(omega, m2, k2, c2, m1, k1, c1):

	X1 = (omega**2 * m1)/(k1 + 1j * omega * c1)

	X2 = (omega**2 * m2)/(k2 + 1j * omega * c2)

	X12 = (omega**2 * m1)/(k2 + 1j * omega * c2)

	X = ((1 - X1)*(1 - X2) - X12)

	VT = np.absolute(1/X)

	return VT



def keff(m1, m2, R, D, mu):
	G = (m1 + m2) * 9.81
	return G/R + mu * G/D


def Xieff(R, D, mu):
	return (2/np.pi) * (mu*R/(D+mu*R))


# Konstanten und Definitionen
m1  = 2486.7
Xi1 = 0.05

m2  = 1619.5
R   = 1.777
D   = 0.325

ag = 3.924
S  = 1
TB = 0.4
TC = 1.6
TD = 2

_T = np.arange(0.001, 4.01, 0.1)
_mu = [0.02, 0.04, 0.10, 0.15, 0.20, 0.25]
AWS = []
AWS_rayleigh = []
for T in _T:
	AWS.append(Se(1, ag, S, TB, TC, TD, T))

for mu in _mu:	
	_AWS = []
	for T in _T:
		k1 = m1 * (2 * np.pi / T)**2
		k2 = keff(m1, m2, R, D, mu)
		Xi2 = Xieff(R, D, mu)
		Sa = vereinfacht_rayleigh(Xi1, k1, m1, Xi2, k2, m2)
		_AWS.append(Sa)
	AWS_rayleigh.append([mu, _AWS])
	

plt.plot(_T, AWS, label='AWS', linestyle='-')

for AWS in AWS_rayleigh:
	plt.plot(_T, AWS[1], label=r'Isolationsspektrum ($\mu = $' + str(AWS[0]) + ')', linestyle='-')


plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=2, frameon=False)

plt.xlabel('T [$s$]')
plt.ylabel('$S_a$ $[m/s^{2}]$')

plt.grid(True)
plt.tight_layout()
#plt.axes().set_aspect(0.7)
plt.margins(x=0, y=0.1)
plt.gca().set_ylim(0,11)
plt.gca().set_xlim(0,4.01)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/Variation.png", dpi=300)
plt.clf()
