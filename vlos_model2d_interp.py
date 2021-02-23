import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pylab as plt
np.warnings.filterwarnings('ignore')
import matplotlib.pylab as plt




 
def Rings(xy_mesh,pa,inc,x0,y0):
	(x,y) = xy_mesh

	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))



	R= np.sqrt(X**2+(Y/np.cos(inc))**2)

	return R




def pix(R,shape,R_min):

	[ny,nx] = shape
	mask = R < R_min
	indices = np.indices((ny,nx))

	pix_y=  indices[0][mask]
	pix_x=  indices[1][mask]


	return pix_x, pix_y




def vlos_interp(xy_mesh,R_depro,Vrot,Vr2,pa,inc,x0,y0,Vsys,Vt2,phi_b,vmode,shape,pixel_scale):



	[ny,nx] = shape
	model = np.zeros((ny,nx))



	pa,inc,phi_b =(pa)*np.pi/180,inc*np.pi/180,phi_b*np.pi/180



	Vrot, Vr2, Vt2 = np.asarray(Vrot),np.asarray(Vr2),np.asarray(Vt2)

	R_depro = np.asarray([])
	for m,n in xy_mesh:
		r = Rings((m,n),pa,inc,x0,y0)
		R_depro = np.append(R_depro, np.median(r))

	min_R = np.nanmin(R_depro)
	#R_depro = np.asarray(R_depro)/pixel_scale
	R_depro = np.insert(R_depro, 0, 0)

	# We need to include R = 0, and V = 0

	if vmode == "circular":

		Vrot_list = np.copy(Vrot)
		Vrot_list = np.insert(Vrot_list, 0, 0)

	if vmode == "radial":
		Vrot_list = np.copy(Vrot)
		Vr2_list = np.copy(Vr2)
		Vrot_list = np.insert(Vrot_list, 0, 0)
		Vr2_list = np.insert(Vr2_list, 0, 0)

	if vmode == "bisymmetric":
		Vrot_list = np.copy(Vrot)
		Vr2_list = np.copy(Vr2)
		Vt2_list = np.copy(Vt2)

		Vrot_list = np.insert(Vrot_list, 0, 0)
		Vr2_list = np.insert(Vr2_list, 0, 0)
		Vt2_list = np.insert(Vt2_list, 0, 0)

	#plt.plot(xy_mesh[0],xy_mesh[1],"ko")
	#plt.show()




	def interp_mod(xy_mesh):
		#global min_R
		#min_R = 2e5

		for x,y in xy_mesh:

			k = 0

			R  = Rings((x,y),pa,inc,x0,y0)

			cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
			sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)
			theta = np.arctan(sin_tetha/cos_tetha)


			if vmode == "circular":
				f0 = interp1d(R_depro, Vrot_list, fill_value = "extrapolate")
				#if np.nanmin(R) < min_R:
				#	min_R = np.nanmin(R)

				Vrot_new = f0(R)
				vlos = np.sin(inc)*Vrot_new*cos_tetha


				for i,j in zip(x,y):
					model[j][i] = vlos[k]
					k = k+1


			if vmode == "radial":

				f0 = interp1d(R_depro, Vrot_list, fill_value = "extrapolate")
				f1 = interp1d(R_depro, Vr2_list, fill_value = "extrapolate")
				Vrot_new = f0(R)
				Vr2_new = f1(R)

				vlos = np.sin(inc)*(Vrot_new*cos_tetha + Vr2_new*sin_tetha)

				for i,j in zip(x,y):
					model[j][i] = vlos[k]
					k = k+1


			if vmode == "bisymmetric":
				f0 = interp1d(R_depro, Vrot_list, fill_value = "extrapolate")
				f1 = interp1d(R_depro, Vr2_list, fill_value = "extrapolate")
				f2 = interp1d(R_depro, Vt2_list, fill_value = "extrapolate")
				Vrot_new = f0(R)
				Vr2_new = f1(R)
				Vt2_new = f2(R)



				m = 2
				theta_b = theta - phi_b
				vlos = np.sin(inc)*((Vrot_new*cos_tetha) - Vt2_new*np.cos(m*theta_b)*cos_tetha - Vr2_new*np.sin(m*theta_b)*sin_tetha)


				for i,j in zip(x,y):
					model[j][i] = vlos[k]
					k = k+1



	interp_mod(xy_mesh)

	#print("out = ", min_R)
	X = np.arange(0, nx, 1)
	Y = np.arange(0, ny, 1)
	XY_mesh = np.meshgrid(X,Y,sparse=True)
	R  = Rings(XY_mesh,pa,inc,x0,y0)
	r = pix(R,shape,min_R)

	interp_mod([[r[0],r[1]]])


	model[model == 0] = np.nan

	return model


