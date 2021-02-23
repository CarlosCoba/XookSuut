import numpy as np
import matplotlib.pylab as plt
import lmfit
import sys
import matplotlib
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
from astropy.stats import sigma_clip
from scipy.interpolate import interp1d
from matplotlib.colors import Normalize
#from numpy.linalg import LinAlgError

from vlos_model2d import vlos
from vlos_model2d_interp import vlos_interp
from poly import legendre



from model_params import M_radial 
from model_params import Rings 
from pixel_params import pixels
import fit_params
from fit_params import fit
from phi_bar_sky import pa_bar_sky


def bisym_mod(vel, evel, guess0, vary, n_it, rstart, rfinal, ring_space, frac_pixel, r_back, delta, pixel_scale, r_bar_max,errors, config, e_ISM):
#def bisym_mod(vel, evel, guess0, vary, n_it=5,rstart = 2, rfinal = 55, ring_space = 2, frac_pixel = 0.7, r_back = 20, delta = 1, pixel_scale = 0.2, r_bar_max = 20):

		vrot0,vr20,pa0,inc0,x0,y0,vsys0,vtan,theta_b = guess0
		vmode = "bisymmetric"
		[ny,nx] = vel.shape
		shape = [ny,nx]


		"""

		 					BYSIMETRIC MODEL


		"""


		chisq_global = 1e4
		PA, INC, XC,YC,VSYS = 0,0,0,0,0
		Vrot, Vrad, Vsys,Vtan = [],[],[],[]
		R = 0

		if r_back < ring_space:
			r_back = ring_space



		for jj in np.arange(0,r_back,ring_space):
			for it in range(n_it):
				guess = [vrot0,vr20,pa0,inc0,x0,y0,vsys0,0,theta_b]

				#theta_b = np.arctan(np.tan(theta_b*np.pi/180-pa0*np.pi/180)/np.cos(inc0*np.pi/180))*180/np.pi

				xi_sq_array = np.asarray([])
				N = np.asarray([])


				vrot_model,vrad_model,vtan_model = np.asarray([]),np.asarray([]),np.asarray([])
				los_vel = np.array([])
				e_los_vel = np.array([])
				x_pos = np.array([])
				y_pos = np.array([])

				los_vel = []
				e_los_vel = []
				x_pos = []
				y_pos = []
				xy_pos = []
				r = []


				for j in np.arange(rstart,rfinal-jj,ring_space):

					XY_mesh, vel_val, e_vel, f_pixel = pixels(vel,evel,guess,ringpos = j, delta=delta,ring = "ARCSEC",pixel_scale=pixel_scale)
					npixels = len(XY_mesh[0])


					if j < 10: 
						f_pixel = 1


					if f_pixel > frac_pixel:


						if j < r_bar_max:

							# Create bisymetric model
							try:
								w_sys,w_rot,w_rad,w_tan,vrot,vrad,vsys,vtan = M_radial(XY_mesh,0,0,pa0,inc0,x0,y0,vsys0,0,theta_b,vel_val-vsys0,e_vel,vmode)
								vrot_model,vrad_model,vtan_model = np.append(vrot_model,vrot),np.append(vrad_model,vrad),np.append(vtan_model,vtan)
							except(np.linalg.LinAlgError):
								vrot_model,vrad_model,vtan_model = np.append(vrot_model,100),np.append(vrad_model,10),np.append(vtan_model,10)							

						else:

							# Create ciruclar model
							w_rot,w_rad,vrot,vrad = M_radial(XY_mesh,0,0,pa0,inc0,x0,y0,vsys0,0,theta_b,vel_val-vsys0,e_vel,vmode= "circular")
							vrot_model,vrad_model,vtan_model = np.append(vrot_model,vrot),np.append(vrad_model,0),np.append(vtan_model,0)

						r.append(j)
						los_vel.append(vel_val)
						e_los_vel.append(e_vel)
						x_pos.append(XY_mesh[0])
						y_pos.append(XY_mesh[1])
						xy_pos.append(XY_mesh)





				guess = [vrot_model+1,vrad_model+1,pa0,inc0,x0,y0,vsys0, vtan_model+1,theta_b]

				vrot , vrad, vsys0,  pa0, inc0 , x0, y0, vtan, theta_b, xi_sq, n_data = fit(shape,los_vel,e_los_vel,xy_pos,guess,vary,vmode,config,fit_method = "Powell", r_bar_max = r_bar_max, e_ISM =  e_ISM)



				if xi_sq < chisq_global:

					PA, INC, XC,YC,VSYS,THETA = pa0, inc0, x0, y0,vsys0,theta_b
					LOS_VEL,eLOS_VEL = los_vel,e_los_vel
					Vrot = vrot
					Vrad = vrad
					Vtan = vtan
					R = r
					XY_PIX_POS = xy_pos
					chisq_global = xi_sq
					#print("chisq_global1 = ", chisq_global)

				if it == n_it -1 :

					vrot_poly = legendre(R,Vrot)
					vrad_poly = legendre(R,Vrad)
					vtan_poly = legendre(R,Vtan)
					guess_end = [vrot_poly,vrad_poly,PA,INC,XC,YC,VSYS,vtan_poly,THETA]
					guess_end = [vrot_poly,Vrad,PA,INC,XC,YC,VSYS,Vtan,THETA]
					vary_end = [True,True,False,False,False,False,False,True,False]
					#vary_end = [True,True,True,True,True,True,True,True,True]

					Vrot , Vrad, vsys0,  pa0, inc0 , x0, y0, Vtan, theta_b, xi_sq, n_data = fit(shape,LOS_VEL,eLOS_VEL,XY_PIX_POS,guess_end,vary_end,vmode,"",fit_method = "Powell", r_bar_max = r_bar_max, e_ISM =  e_ISM)

					#print("chisq_global2=", xi_sq)

					MODEL_not_interp = np.zeros((ny,nx)) 
					N = len(LOS_VEL)
					for k in range(N):
						for mm,nn in zip(XY_PIX_POS[k][0],XY_PIX_POS[k][1]): 
							MODEL_not_interp[nn][mm] = vlos(mm,nn,Vrot[k],Vrad[k],PA,INC,XC,YC,VSYS,Vtan[k],THETA,vmode) - VSYS

					MODEL_interp = vlos_interp(XY_PIX_POS,R,Vrot,Vrad,PA,INC,XC,YC,VSYS,Vtan,THETA,vmode,shape,pixel_scale)

					# For error estimation:
					LoS,eLoS,XY_pixels =  los_vel,e_los_vel,xy_pos
					best = guess





		MODEL_not_interp[MODEL_not_interp == 0] = np.nan
		MODELS = [MODEL_interp,MODEL_not_interp]


		R = np.asarray(R)
		Vtan = np.asarray(Vtan)
		Vrad = np.asarray(Vrad)
		Vrot = np.asarray(Vrot)
		PA_bar_major = pa_bar_sky(PA,INC,THETA)
		PA_bar_minor = pa_bar_sky(PA,INC,THETA-90)

		if errors == 1:
			from mcmc import fit_mcmc
			res_mcmc =  fit_mcmc(shape, LOS_VEL, eLOS_VEL, XY_PIX_POS, guess_end, vary_end, vmode, "", fit_method = "emcee", e_ISM = e_ISM)
		else:
			res_mcmc = Vrot*0,[Vrot*0,Vrot*0],Vrot*0,[Vrot*0,Vrot*0],0,[0,0],0,[0,0],0,[0,0],0,[0,0],0,[0,0],0,[0,0],Vrot*0,[Vrot*0,Vrot*0]

		return PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODELS,PA_bar_major,PA_bar_minor,chisq_global,res_mcmc
