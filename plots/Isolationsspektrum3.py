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
def Se(ag, S, TB, TC, TD, T):

	if T <= TB:
		Se = ag * S * (1 + T/TB * (2.5 -1))

	elif T <= TC:
		Se = ag * S * 2.5

	elif T <= TD:
		Se = ag * S * 2.5 * TC/T

	else:
		Se = ag * S * 2.5 * (TC*TD)/(T*T)

	return Se

# Übertagung von Bodenanregung zu Masse 1
def VT_1_over_3(omega, Xi1, m2, k2, c2):

	# Parameter des EMS aus dem AWS
	m1 = 1
	k1 = omega**2
	c1 = Xi1 * 2 * np.sqrt(k1)


	X1 = (omega**2 * m1)/(k1 + 1j * omega * c1)

	X2 = (omega**2 * m2)/(k2 + 1j * omega * c2)

	X12 = (omega**2 * m1)/(k2 + 1j * omega * c2)

	X = ((1 - X1)*(1 - X2) - X12)

	VT = np.absolute(1/X)

	return VT


# Übertagung von Bodenanregung zu Masse 2
def VT_2_over_3(omega, Xi1, m2, k2, c2):

	# Parameter des EMS aus dem AWS
	m1 = 1
	k1 = omega**2
	c1 = Xi1 * 2 * np.sqrt(k1)


	X1 = (k1 + 1j * omega * c1) / (k1 + 1j * omega * c1 - omega**2 * m1)

	X2 = (omega**2 * m2)/(k2 + 1j * omega * c2)

	X3 = (k1 + 1j * omega * c1)/(k2 + 1j * omega * c2)

	X = 1 - X1*X3 - X2 - X3

	VT = np.absolute(1/X)

	return VT


# Beteiligungsfaktor
def Phi(m1, k1, m2, k2):

	A = (k2 + k1) * m1 + k1 * m2
	B = ((k2 + k1) * m1 + k1 * m2)**2 - 4 * m2 * m1 * k2 * k1
	C = 2 * m2 * m1

	omega1 = np.sqrt((A - np.sqrt(B))/C)
	phi = (k2 + k1 - m2 * omega1**2)/k1

	return 1/phi



# Konstanten und Definitionen
m1  = 1
k1  = 1000
Xi1 = 0.05

k2  = 200
Xi2 = 0.05

ag = 0.4
S  = 0.75
TB = 0.1
TC = 0.5
TD = 2

_mu = np.arange(1, 100, 1)
F_Kelly = []
F_Iso   = []
F_Iso2  = []

for mu in _mu:

	m2 = m1 * mu
	c2 = Xi2 * 2 * np.sqrt(k2 * m2)

	# Kelly
	omega = np.sqrt(k2 / (m1 + m2))
	T = (2 * np.pi ) / omega

	F_Kelly.append(Se(ag, S, TB, TC, TD, T) * (m1 + m2) * Phi(m1, k1, m2, k2))


	# Isolationsspektrum 1
	omega = np.sqrt(k1 / m1)
	T = (2 * np.pi ) / omega

	F_Iso.append(Se(ag, S, TB, TC, TD, T) * VT_1_over_3(omega, Xi1, m2, k2, c2) * m1 )


	# Isolationsspektrum 2
	omega = np.sqrt(k1 / m1)
	T = (2 * np.pi ) / omega

	F_Iso2.append(Se(ag, S, TB, TC, TD, T) * VT_2_over_3(omega, Xi1, m2, k2, c2) * m1 )





	


plt.plot(_mu, F_Kelly, label='Kelly', linestyle='-')
plt.plot(_mu, F_Iso, label='Isoliert S1/S3', linestyle=':')
plt.plot(_mu, F_Iso2, label='Isoliert S2/S3', linestyle='--')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=3, frameon=False)

plt.xlabel('m2/m1')

plt.grid(True)
plt.tight_layout()
#plt.axes().set_aspect(0.7)
plt.margins(x=0, y=0)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/Isolation.png", dpi=300)
plt.clf()
