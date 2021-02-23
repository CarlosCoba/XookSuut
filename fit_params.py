import numpy as np
import numpy as np
import matplotlib.pylab as plt
import scipy
import sys
import lmfit
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
from matplotlib.gridspec import GridSpec
#import sys

from recompute_chi import result
from read_config import config_file
 
def Rings(xy_mesh,pa,inc,x0,y0):
	(x,y) = xy_mesh

	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))



	R= np.sqrt(X**2+(Y/np.cos(inc))**2)

	return R





def RADIAL_MODEL(xy_mesh,Vrot,Vr2,pa,inc,x0,y0,Vsys):
	(x,y) = xy_mesh
	pa,inc=(pa)*np.pi/180,inc*np.pi/180
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)
	vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha + Vr2*sin_tetha)
	return np.ravel(vlos)


def CIRC_MODEL(xy_mesh,Vrot,pa,inc,x0,y0,Vsys):
	(x,y) = xy_mesh
	pa,inc=(pa)*np.pi/180,inc*np.pi/180
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)
	vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha)
	return np.ravel(vlos)




def RADIAL_MODEL2(w1,Vsys,w2,Vrot,w3,Vrad):
	vlos = w1*Vsys + w2*Vrot + w3*Vrad
	return np.ravel(vlos)


def BISYM_MODEL(xy_mesh,Vrot,Vrad,pa,inc,x0,y0,Vsys,Vtan,phi_b):

	(x,y) = xy_mesh
	pa,inc,phi_b=(pa)*np.pi/180,inc*np.pi/180,phi_b*np.pi/180
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)

	m = 2
	theta = np.arctan(sin_tetha/cos_tetha)
	#phi_b = np.arctan(np.tan(phi_b-pa)/np.cos(inc))

	phi_b = phi_b
	theta_b = theta - phi_b
	

	vlos = Vsys+np.sin(inc)*((Vrot*cos_tetha) - Vtan*np.cos(m*theta_b)*cos_tetha - Vrad*np.sin(m*theta_b)*sin_tetha)

	return np.ravel(vlos)



def polynomial(x,a0,a1,a2,a3,a4,a5):
	x = np.asarray(x)
	y = a0 + a1*x +a2*x**2 +a3*x**3 + a4*x**4 + a5*x**5
	return y

def linear(x,a0,a1):
	x = np.asarray(x)
	y = a0 + a1*x 
	return y

def fit_polynomial(x,dato):
	x = np.asarray(x)
	dato = np.asarray(dato)


	def residual_line(pars,x,data=None):
		parvals = pars.valuesdict()
		a0 = parvals['a0']
		a1 = parvals['a1']
		a2 = parvals['a2']
		a3 = parvals['a3']
		a4 = parvals['a4']
		a5 = parvals['a5']
		model = polynomial(x,a0,a1,a2,a3,a4,a5)
		objective = model - data
		return objective**2


	fit_param = Parameters()
	fit_param.add('a0', value=0)
	fit_param.add('a1', value=0)
	fit_param.add('a2', value=0)
	fit_param.add('a3', value=0)
	fit_param.add('a4', value=0)
	fit_param.add('a5', value=0)

	out = minimize(residual_line, fit_param, args=(x,), kws={'data': dato},method='Powell', nan_policy = "omit")
	best = out.params
	a0, a1, a2, a3, a4, a5 = best["a0"].value,best["a1"].value, best["a2"].value, best["a3"].value, best["a4"].value, best["a5"].value
	return a0, a1, a2, a3, a4, a5


def fit_linear(x,dato):
	x = np.asarray(x)
	dato = np.asarray(dato)


	def residual_line(pars,x,data=None):
		parvals = pars.valuesdict()
		a0 = parvals['a0']
		a1 = parvals['a1']
		model = linear(x,a0,a1)
		objective = model - data
		return objective**2


	fit_param = Parameters()
	fit_param.add('a0', value=1e3)
	fit_param.add('a1', value=0)

	out = minimize(residual_line, fit_param, args=(x,), kws={'data': dato},method='Powell')
	best = out.params
	a0, a1 = best["a0"].value,best["a1"].value
	return a0, a1

def fit(shape,vel_val,e_vel,xy_mesh,guess,vary,vmode,config,fit_method = "leastsq",e_ISM = 5, r_bar_max = 15):
	"""
	vary = [Vrot,Vrad,PA,INC,XC,YC,VSYS,theta,Vtan]
	"""
	vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0 = guess
	constant_params = [pa0,inc0,X0,Y0,vsys0,theta0]


	[ny,nx] = shape
	if vmode == "circular":



#####
		def vmodel_dataset(pars, xy_mesh, i):

			parvals = pars.valuesdict()
			pa = parvals['pa']
			inc = parvals['inc']
			#Vsys = parvals['Vsys_%i'% (i+1)]
			Vsys = parvals['Vsys']
			Vrot = parvals['Vrot_%i'% (i+1)]
			x0,y0 = parvals['x0'],parvals['y0']

			model = CIRC_MODEL(xy_mesh,Vrot,pa,inc,x0,y0,Vsys)

			return model

		def residual(pars, xy_mesh, data= None, weight = e_vel):
			"""Calculate total residual for fits of VMODELS to several data sets."""

			ndata = len(data)
			resid = np.array([])


			#parvals = pars.valuesdict()
			#e_ISM = parvals['e_ISM']



			# make residual per data set
			for i in range(ndata):

				sigma = np.sqrt(weight[i]**2 + e_ISM**2)
				Resid = (data[i] - vmodel_dataset(pars, xy_mesh[i], i))/sigma
				resid = np.append(resid,Resid)
				
			# now flatten this to a 1D array, as minimize() needs

			return resid.flatten()




    
#####



		fit_params = Parameters()

		for iy, y in enumerate(vel_val):
				fit_params.add('Vrot_%i' % (iy+1),value=vrot0[iy], vary = vary[0], min = 0, max = 400)

		if config == "":

			fit_params.add('Vsys', value=vsys0, vary = vary[6])
			fit_params.add('pa', value=pa0, vary = vary[2], min = 0, max = 360)
			fit_params.add('inc', value=inc0, vary = vary[3], min = 0, max = 90)
			fit_params.add('x0', value=X0, vary = vary[4],  min = 0, max = nx)
			fit_params.add('y0', value=Y0, vary = vary[5], min = 0, max = ny)


		else:
			k = 0
			for res,const in zip(config_file(config),constant_params):
				#param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), float(res["val"]), float(res["min"]), float(res["max"]) 
				param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), eval(res["val"]), eval(res["min"]), eval(res["max"]) 
				fit_params.add(param, value = constant_params[k], vary = fit, min = vmin, max = vmax)
				k = k+1


		out = minimize(residual, fit_params, args=(xy_mesh,), kws={'data': vel_val},method = fit_method, nan_policy = "omit")
		best = out.params

		N = len(vel_val)
		Vsys, Vrad, Vrot = [],[],[]
		for iy in range(1,N+1):
				Vrot.append(best["Vrot_%s"%iy].value)



		pa, inc = best["pa"].value,best["inc"].value
		x0,y0  = best["x0"].value,best["y0"].value
		Vsys = best["Vsys"].value

		N_free = out.nfree
		red_chi = out.redchi



		#Vlist = [Vrot,0,0]
		#R_v = sigma_R
		#Chisq = result(vmode,vel_val,xy_mesh,e_vel,Vlist,R_v,pa,inc,x0,y0,Vsys,N_free)
		#print(Chisq/N_free, red_chi)
		#red_chi = Chisq


		if inc > 85: inc = 80
		if inc < 20: inc = 30

		if pa > 359: pa = pa -10

		return Vrot, Vsys, pa, inc , abs(x0), abs(y0) ,red_chi, N_free

	if vmode == "radial":

		vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0 = guess


		#XC,YC,V0,Pos_Ang,inclination,theta=[],[],[],[],[],[]
		#rings = []



#####
		def vmodel_dataset(pars, xy_mesh, i):

			parvals = pars.valuesdict()
			pa = parvals['pa']
			inc = parvals['inc']
			#Vsys = parvals['Vsys_%i'% (i+1)]
			Vsys = parvals['Vsys']
			Vrot = parvals['Vrot_%i'% (i+1)]
			Vr2 = parvals['Vr2_%i'% (i+1)]
			x0,y0 = parvals['x0'],parvals['y0']

			model = RADIAL_MODEL(xy_mesh,Vrot,Vr2,pa,inc,x0,y0,Vsys)
			return model

		def residual(pars, xy_mesh, data= None, weight = e_vel):
			"""Calculate total residual for fits of VMODELS to several data sets."""

			ndata = len(data)
			resid = np.array([])


			#parvals = pars.valuesdict()
			#e_ISM = parvals['e_ISM']


			# make residual per data set
			for i in range(ndata):
				sigma = np.sqrt(weight[i]**2 + e_ISM**2)
				Resid = (data[i] - vmodel_dataset(pars, xy_mesh[i], i))/sigma
				#print("Resid = ", Resid)
				resid = np.append(resid,Resid)
				
			# now flatten this to a 1D array, as minimize() needs

			return resid.flatten()




    
#####




		fit_params = Parameters()

		for iy, y in enumerate(vel_val):

				fit_params.add('Vrot_%i' % (iy+1),value=vrot0[iy], vary = vary[0], min = 5, max = 400)
				fit_params.add('Vr2_%i' % (iy+1), value=vr20[iy], vary = vary[1], min = -200, max = 200)

		if config == "":

		
			fit_params.add('Vsys', value=vsys0, vary = vary[6])
			fit_params.add('pa', value=pa0, vary = vary[2], min = 0, max = 360)
			fit_params.add('inc', value=inc0, vary = vary[3], min = 0, max = 90)
			fit_params.add('x0', value=X0, vary = vary[4], min = 0, max = nx)
			fit_params.add('y0', value=Y0, vary = vary[5], min = 0, max = ny)

		else:
			k = 0
			for res,const in zip(config_file(config),constant_params):
				#param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), float(res["val"]), float(res["min"]), float(res["max"]) 
				param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), eval(res["val"]), eval(res["min"]), eval(res["max"]) 
				fit_params.add(param, value = constant_params[k], vary = fit, min = vmin, max = vmax)
				k = k+1



		out = minimize(residual, fit_params, args=(xy_mesh,), kws={'data': np.ravel(vel_val)},method = fit_method, nan_policy = "omit")
		best = out.params

		N = len(vel_val)
		Vsys, Vrad, Vrot = [],[],[]
		for iy in range(1,N+1):
				Vrot.append(best["Vrot_%s"%iy].value)
				Vrad.append(best["Vr2_%s"%iy].value)



		pa, inc = best["pa"].value,best["inc"].value
		x0,y0  = best["x0"].value,best["y0"].value
		Vsys = best["Vsys"].value

		N_free = out.nfree
		red_chi = out.redchi

		if inc > 85: inc = 80
		if inc < 20: inc = 30


		if pa > 359: pa = pa -10

		return Vrot , Vrad, Vsys, pa, inc , abs(x0), abs(y0) ,red_chi, N_free


	if vmode == "bisymmetric":

		vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0 = guess



		def vmodel_dataset(pars, xy_mesh, i, r_bar_max):

			parvals = pars.valuesdict()
			pa = parvals['pa']
			inc = parvals['inc']
			Vtan = parvals['Vtan_%i'% (i+1)]
			Vsys = parvals['Vsys']
			Vrot = parvals['Vrot_%i'% (i+1)]
			Vrad = parvals['Vrad_%i'% (i+1)]
			x0,y0 = parvals['x0'],parvals['y0']
			theta = parvals['phi_b']


			R = Rings(xy_mesh,pa,inc,x0,y0)
			mask_R = R < r_bar_max

			model = BISYM_MODEL(xy_mesh,Vrot,Vrad,pa,inc,x0,y0,Vsys,Vtan,theta)
			
			return model

		def residual_bisym(pars, xy_mesh, data= None, weight = e_vel, r_bar_max=10):
			"""Calculate total residual for fits of VMODELS to several data sets."""

			ndata = len(data)
			resid = np.array([])

			#parvals = pars.valuesdict()
			#e_ISM = parvals['e_ISM']


			# make residual per data set
			for i in range(ndata):

				sigma = np.sqrt(weight[i]**2 + e_ISM**2)
				Resid = (data[i] - vmodel_dataset(pars, xy_mesh[i], i, r_bar_max ))/sigma
				resid = np.append(resid,Resid)
				
			# now flatten this to a 1D array, as minimize() needs

			return resid.flatten()



		N = len(vel_val)
		fit_params = Parameters()


		for iy, y in enumerate(vel_val):
				if vr20[iy]-1 == 0 and vt20[iy]-1 ==0:
					vary_vrad = False
					vary_vtan = False


				else:

					vary_vrad = vary[1]
					vary_vtan = vary[7]
					#vary_theta = vary[8]


				fit_params.add('Vrot_%i' % (iy+1),value=vrot0[iy], vary = vary[0], min = 0, max = 400)
				fit_params.add('Vrad_%i' % (iy+1), value=vr20[iy], vary = vary_vrad,  min = -200, max = 200)
				fit_params.add('Vtan_%i' % (iy+1), value=vt20[iy], vary = vary_vtan, min = -200, max = 200)




		if config == "":
			fit_params.add('Vsys', value=vsys0, vary = vary[6])
			fit_params.add('pa', value=pa0, vary = vary[2], min = 0, max = 360)
			fit_params.add('inc', value=inc0, vary = vary[3], min = 0, max = 90)
			fit_params.add('x0', value=X0, vary = vary[4], min = 0, max = nx)
			fit_params.add('y0', value=Y0, vary = vary[5], min = 0, max = nx)
			fit_params.add('phi_b', value=theta0, vary = vary[7], min = 0 , max = 180)
	
		else:
			k = 0
			for res,const in zip(config_file(config),constant_params):
				#param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), float(res["val"]), float(res["min"]), float(res["max"]) 
				param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), eval(res["val"]), eval(res["min"]), eval(res["max"]) 
				fit_params.add(param, value = constant_params[k], vary = fit, min = vmin, max = vmax)
				k = k+1


		out = minimize(residual_bisym, fit_params, args=(xy_mesh,), kws={'data': vel_val, "r_bar_max": r_bar_max},method = fit_method, nan_policy = "omit")
		best = out.params


		N = len(vel_val)

		Vtan, Vrad, Vrot = [],[],[]
		for iy in range(1,N+1):
				Vtan.append(best["Vtan_%s"%iy].value)
				Vrot.append(best["Vrot_%s"%iy].value)
				Vrad.append(best["Vrad_%s"%iy].value)



		pa, inc = best["pa"].value,best["inc"].value
		x0,y0  = best["x0"].value,best["y0"].value
		Vsys = best["Vsys"].value
		theta = best["phi_b"].value

		if theta > 179 or theta < -179: theta = 45
		if theta > 0 and theta < 1: theta = 45

		N_free = out.nfree
		red_chi = out.redchi



		#Vlist = [Vrot,Vrad,Vtan]
		#R_v = sigma_R
		#Chisq = result(vmode,vel_val,xy_mesh,e_vel,Vlist,R_v,pa,inc,x0,y0,Vsys,N_free, theta)
		#print(Chisq/N_free, red_chi,out.chisqr/N_free)
		#red_chi = Chisq
		#print()




		if inc > 85: inc = 80
		if inc < 20: inc = 30

		return Vrot , Vrad, Vsys,  pa, inc , x0, y0, Vtan, theta, red_chi,N_free







