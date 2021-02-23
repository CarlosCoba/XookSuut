import numpy as np
import matplotlib.pylab as plt
import sys
import matplotlib
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

from scipy import special
from scipy.special import gamma
import random



def xi_sq_pdf(x,k=3):
	xi = (1/2.)**(k/2.)*(np.exp(-x/2.))*x**(k/2.-1)/gamma(k/2.)
	return xi

 
def Rings(nx,ny,pa,inc,x0,y0):

	X = np.arange(0, nx, 1)
	Y = np.arange(0, ny, 1)
	XY_mesh = np.meshgrid(X,Y,sparse=True)
	(x,y) = XY_mesh

	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))



	R= np.sqrt(X**2+(Y/np.cos(inc))**2)

	return R


def arctan(r_sky,v_max,R_turn):
	return (2/np.pi)*v_max*np.arctan(r_sky/R_turn)


def vrot(r,vmax,Rturn):
	x = r/Rturn
	return (2/np.pi)*vmax*np.arctan(x)


# Bertola 1991
def bertola(r_sky,v_max,kk,Gamma):
	v = v_max*r_sky/(r_sky**2 + kk**2)**(Gamma/2.)
	return v



#From Haeun Chung, to model the rising RC
def tanh_linear(r_sky,v_max,R1,R2):
		R_turn = R1
		# R2 can be negatica in case of descending RCs
		v=v_max*(np.tanh(r_sky/R_turn) + r_sky*R2 )
		return v


def save_fits(array,name, header = [], value = []):


		hdu = fits.ImageHDU()
		hdu.data = array

		if len(header) != 0:
			for hd,vl in zip(header,value):

				hdu.header[hd] = vl


		hdu.header['FILE0'] = 'VLOS_model'
		hdu.header['UNITS0'] = "km/s"
		hdu.writeto("./models_califa/%s.fits"%name,overwrite=True)




def disk_fit_models(nx,ny,v_max,kk,Gamma,pa,inc,x0,y0,Vsys,phi_bar,save = False,name = False,**kwds):

	print(pa,inc,phi_bar,v_max,kk,Gamma)
	#print(pa,inc,phi_bar,v_max,R1, R2)
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	phi_bar = phi_bar*np.pi/180
	#R  = Rings(nx,ny,pa,inc,x0,y0)

	X = np.arange(0, nx, 1)
	Y = np.arange(0, ny, 1)
	XY_mesh = np.meshgrid(X,Y,sparse=True)
	(x,y) = XY_mesh

	X = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))
	Y = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))



	R= np.sqrt(X**2+(Y/np.cos(inc))**2)




	Vrot = arctan(R,170,7)
	Vrot = bertola(R,v_max,kk,Gamma)
	#Vrot = tanh_linear(R,v_max,R1,R2)
	#Vrot = bertola(R,200,3,1.4)
	Vt2 = xi_sq_pdf(R,8.6)*300
	Vr2 = xi_sq_pdf(R,7.3)*250



	plt.plot(R,Vt2,"dodgerblue")
	plt.plot(R,Vr2,"orange")
	plt.plot(R,Vrot,"k")
	plt.xlim(0,50)
	plt.show()


	#theta_b = np.arctan(np.tan(phi_bar - PA)/np.cos(inc))
	

	cos_tetha = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))/R
	sin_tetha = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))/(np.cos(inc)*R)
	theta = np.arctan(sin_tetha/cos_tetha)
	theta_b = theta - phi_bar

	vlos = Vsys+np.sin(inc)*((Vrot*cos_tetha)-Vt2*np.cos(2*theta_b)*cos_tetha - Vr2*np.sin(2*theta_b)*sin_tetha)
	#vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha)
	#vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha + Vr2*sin_tetha)

	e_ISM, e_D = 2,2
	mu,sigma = 0,np.sqrt(e_ISM**2 + e_D**2)
	noise =  np.random.normal(mu, sigma, size = (ny,nx))


	fwhm = 2.5
	std = fwhm/2.354
	kernel = Gaussian2DKernel(x_stddev=std)
	VLOS = convolve(vlos, kernel)
	VLOS = vlos + noise

	if save == True:
		save_fits(VLOS,name,**kwds)


	plt.imshow(VLOS,origin = "l")
	plt.colorbar()
	plt.show()
	return VLOS


for n in range(10):
	nx,ny = 77,72
	inc = random.randint(20,85)
	PA = random.randint(0,360)
	phi_bar = random.randint(0,180)
	v_max,kk,Gamma = random.randint(100,320),random.uniform(0,60),random.uniform(1,1.5)
	#v_max,R1,R2 = random.randint(80,320),random.uniform(0,60),random.uniform(0,10)
	pa,inc,x0,y0,vsys,phi_bar = PA,inc,36.5,32.5,490,phi_bar
	model = disk_fit_models(nx,ny,v_max,kk,Gamma,pa,inc,x0,y0,vsys,phi_bar, save = False,header = ["PA", "INC","XC","YC","VSYS","PHI_BAR"], value = [pa,inc,x0,y0,vsys,phi_bar],name = "test2_bm.califa_%s"%n)
	#model = disk_fit_models(nx,ny,v_max,R1,R2,pa,inc,x0,y0,vsys,phi_bar, save = False,header = ["PA", "INC","XC","YC","VSYS","PHI_BAR"], value = [pa,inc,x0,y0,vsys,phi_bar],name = "test2_bm.califa_%s"%n)


