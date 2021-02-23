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
from pixel_params import pixels
import fit_params
from fit_params import fit




def circ_mod(vel, evel, guess0, vary, n_it, rstart, rfinal, ring_space, frac_pixel, r_back, delta, pixel_scale, r_bar_max, errors, config, e_ISM ):
#def circ_mod(vel, evel, guess0, vary, n_it=5,rstart = 2, rfinal = 55, ring_space = 2, frac_pixel = 0.7, r_back = 20, delta = 1, pixel_scale = 0.2, r_bar_max = 20):

		vrot0,vr20,pa0,inc0,x0,y0,vsys0,vtan,theta_b = guess0
		vmode = "circular"
		[ny,nx] = vel.shape
		shape = [ny,nx]

		"""

		 					CIRCULAR MODEL


		"""

		chisq_global = 1e10
		PA, INC, XC,YC,VSYS = 0,0,0,0,0
		Vrot, Vrad, Vsys,Vtan = [],[],[],[]
		R = 0
		best_xy_pix = []

		if r_back < ring_space:
			r_back = ring_space


		for jj in np.arange(0,r_back,ring_space):
			for it in np.arange(n_it):
				guess = [vrot0,vr20,pa0,inc0,x0,y0,vsys0,0,theta_b]

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


					if f_pixel > frac_pixel:

						# Create model
						try:
							w_rot,w_rad,vrot,vrad = M_radial(XY_mesh,0,0,pa0,inc0,x0,y0,vsys0,0,0,vel_val-vsys0,e_vel,vmode)
							vrot_model,vrad_model = np.append(vrot_model,vrot),np.append(vrad_model,vrad)



						except(np.linalg.LinAlgError):
							vrot_model,vrad_model = np.append(vrot_model,100),np.append(vrad_model,0)



						r.append(j)

						los_vel.append(vel_val)
						e_los_vel.append(e_vel)
						x_pos.append(XY_mesh[0])
						y_pos.append(XY_mesh[1])
						xy_pos.append(XY_mesh)




				guess = [vrot_model+1,vrad_model+1,pa0,inc0,x0,y0,vsys0, vtan_model,theta_b]
				vrot , vsys0,  pa0, inc0, x0, y0, xi_sq, n_data = fit(shape, los_vel, e_los_vel, xy_pos, guess, vary, vmode, config, fit_method = "Powell", e_ISM = e_ISM)

				if xi_sq < chisq_global:

					LOS_VEL = los_vel
					e_LOS_VEL = e_los_vel
					PA, INC, XC,YC,VSYS,THETA = pa0, inc0, x0, y0,vsys0,theta_b
					R = r
					Vrot = vrot
					XY_PIX_POS = xy_pos
					chisq_global = xi_sq

					#print("chisq_global1 = ", chisq_global)


				if it == n_it -1 :
					


					vrot_poly = legendre(R,Vrot)
					vary_end = [True,True,False,False,False,False,False,True,True]
					guess_end = [vrot_poly,vr20,PA,INC,XC,YC,VSYS,0,THETA]
					Vrot , vsys0,  pa0, inc0, x0, y0, xi_sq, n_data = fit(shape, LOS_VEL, e_LOS_VEL, XY_PIX_POS, guess_end, vary_end, vmode, "", fit_method = "Powell", e_ISM = e_ISM)
					#print("chisq_global2 = ", xi_sq)

					MODEL_not_interp = np.zeros((ny,nx)) 
					N = len(LOS_VEL)
					for k in range(N):
						for mm,nn in zip(XY_PIX_POS[k][0],XY_PIX_POS[k][1]): 
							MODEL_not_interp[nn][mm] = vlos(mm,nn,Vrot[k],0,PA,INC,XC,YC,VSYS,0,0,vmode) - VSYS

					MODEL_interp = vlos_interp(XY_PIX_POS,R,Vrot,0,PA,INC,XC,YC,VSYS,0,0,vmode,shape,pixel_scale)

					# For error estimation:
					LoS,eLoS,XY_pixels =  los_vel,e_los_vel,xy_pos
					best = guess



		MODEL_not_interp[MODEL_not_interp == 0] = np.nan
		MODELS = [MODEL_interp,MODEL_not_interp]
		Vrot = np.array(Vrot)
		R = np.array(R)


		if errors == 1:

			from mcmc import fit_mcmc
			res_mcmc =  fit_mcmc(shape, LOS_VEL, e_LOS_VEL, XY_PIX_POS, guess_end, vary_end, vmode, "", fit_method = "emcee", e_ISM = e_ISM)
		else:
			res_mcmc = Vrot*0,[Vrot*0,Vrot*0],Vrot*0,[Vrot*0,Vrot*0],0,[0,0],0,[0,0],0,[0,0],0,[0,0],0,[0,0],0,[0,0],Vrot*0,[Vrot*0,Vrot*0]

		return PA,INC,XC,YC,VSYS,0,R,Vrot,0*Vrot,0*Vrot,MODELS,chisq_global,res_mcmc


