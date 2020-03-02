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


# Kelly
def Kelly(n, m1, k2, m2):

	omega = np.sqrt(k2/(m1+m2))
	T1 = 2 * np.pi / omega

	return Se(n, ag, S, TB, TC, TD, T1)


# Isemann
def vereinfacht(n, k1, m1, k2, m2):

	omega = np.sqrt(k1/m1)

	a = k2/k1
	b = m2/m1

	A = 1 + a + b
	B = (1 + a + b)**2 - 4 * a *b
	C = 2 * b

	omega1 = np.sqrt((A - np.sqrt(B))/(C) * omega**2)
	T1 = 2*np.pi/omega1

	A = 1 + a - b
	B = (1 + a + b)**2 - 4 * a *b

	r = (A + np.sqrt(B))/(2)

	phi = 1/r

	L = (phi * m2 + 1 * m1) / (phi**2 * m2 + 1**2 * m1 )

	return Se(n, ag, S, TB, TC, TD, T1) * phi * L


# Rick
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
	#T1 = 2*np.pi/omega_1d
	T1 = 2*np.pi/omega1

	# Transmissionskoeffizient des Systems

	# Nach Pocanschi
	#a = (2*_omega1*_omega2*(Xi2 * _omega2 - Xi1 * _omega1))/(_omega2**2 - _omega1**2)
	#b = (2*(Xi1 * _omega2 - Xi2 * _omega1))/(_omega2**2 - _omega1**2)
	a = (2*_omega1*_omega2*(Xi2 * _omega2 - Xi1 * _omega1))/(_omega2**2 - _omega1**2)
	b = (2*(Xi1 * _omega2 - Xi2 * _omega1))/(_omega2**2 - _omega1**2)

	#  Nach Isemann
	#b = 2 * Xi2 / (omega1 * omega2)
	#a = omega1 * omega2 * b

	#_c1 = a * _m1 + b * _k1
	#_c2 = a * _m2 + b * _k2
	
	c1 = Xi1 * 2 * np.sqrt(m1*k1)
	c2 = Xi2 * 2 * np.sqrt(m2*k2)

	VT = VT_1_over_3(omega_1d, m2, k2, c2, m1, k1, c1) 

	n = 1/(2*Xi2*10)#np.sqrt(10/5)

	return Se(n, ag, S, TB, TC, TD, T1) * VT / 10.05


def VT_1_over_3(omega, m2, k2, c2, m1, k1, c1):

	X1 = (omega**2 * m1)/(k1 + 1j * omega * c1)

	X2 = (omega**2 * m2)/(k2 + 1j * omega * c2)

	X12 = (omega**2 * m1)/(k2 + 1j * omega * c2)

	X = ((1 - X1)*(1 - X2) - X12)

	VT = np.absolute(1/X)

	return VT


# Übertagung von Bodenanregung zu Masse 2
def VT_2_over_3(omega, m2, k2, c2, m1, k1, c1):

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
m1  = 2486.7
k1  = 98170
Xi1 = 0.05
c1 = Xi1 * 2 * np.sqrt(k1 * m1)

m2  = 1619.5
k2  = 26848.2
Xi2 = 0.14147

ag = 3.924
S  = 1
TB = 0.4
TC = 1.6
TD = 2

_T = np.arange(0.001, 4.01, 0.1)
#_T = [1]
AWS = []
for T in _T:
	AWS.append(Se(1, ag, S, TB, TC, TD, T))
	
AWS_rayleigh = []
for T in _T:
	k1 = m1 * (2 * np.pi / T)**2
	Sa = vereinfacht_rayleigh(Xi1, k1, m1, Xi2, k2, m2)
	AWS_rayleigh.append(Sa)	

AWS_iso = []
for T in _T:
	n = 1/(2*Xi2*10)#np.sqrt(10/(5+Xi2*100))
	k1 = m1 * (2 * np.pi / T)**2
	Sa = vereinfacht(n, k1, m1, k2, m2)
	AWS_iso.append(Sa)

AWS_n = []
for T in _T:
	n = 1/(2*Xi2*10)#np.sqrt(10/(5+Xi2*100))
	AWS_n.append(Se(n, ag, S, TB, TC, TD, T))

plt.plot(_T, AWS, label='AWS', linestyle='-')
plt.plot(_T, AWS_n, label=r'AWS ($\eta = 0.335$)', linestyle='-')
T_iso = 2 * np.pi / np.sqrt(k2/(m2+m1))
#plt.axvline(x=T_iso, label = 'Isolatorperiode', color='black', linestyle='--', linewidth=0.5)

plt.plot(_T, AWS_rayleigh, label='Isolationsspektrum (Transmissibilität)', linestyle=':')
plt.plot(_T, AWS_iso, label='AWS isoliert (vereinfacht)', linestyle='--')

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=2, frameon=False)

plt.xlabel('T [$s$]')
plt.ylabel('$S_a$ $[m/s^{2}]$')

plt.grid(True)
plt.tight_layout()
#plt.axes().set_aspect(0.7)
plt.margins(x=0, y=0.1)
plt.gca().set_ylim(0,11)
plt.gca().set_xlim(0,4.01)
plt.savefig("/home/couchsofa/masterthesis/Masterthesis/images/Isolation_5.png", dpi=300)
plt.clf()