import numpy as np
import matplotlib.pylab as plt


def plot_true(vmode,axis):




	R = np.linspace(0,60,200)
	Vrot = 80*(R/27.8)**0.35
	R1 = R[R < 10*np.pi]
	R1[R1 == 0] = np.nan
	Vt2 = 25*np.sin(R1/10.)
	Vr2 = 0.85*Vt2

	if vmode == "bisymmetric": 

		axis.plot(R1,Vt2,color = "dodgerblue",linestyle='-', alpha = 0.5, linewidth=2)
		axis.plot(R1,Vr2,color = "orange",linestyle='-', alpha = 0.5, linewidth=2)
		axis.plot(R,Vrot,color = "k",linestyle='-', alpha = 0.5, linewidth=2)

	if vmode == "radial": 

		axis.plot(R1,Vr2,color = "orange",linestyle='-', alpha = 0.5, linewidth=2)
		axis.plot(R,Vrot,color = "k",linestyle='-', alpha = 0.5, linewidth=2)

	if vmode == "circular":
		axis.plot(R,Vrot,color = "k",linestyle='-', alpha = 0.5, linewidth=2)

