import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
np.warnings.filterwarnings('ignore')
import os.path
from os import path

from axis import AXIS 
from fit_rotcurve import Rturn_Vmax


from read_config import config_file
from xs import init_models 

from kinematic_centre_vsys import KC
from CBAR import colorbar as cb

from isolated_pixels import filter_isolated_cells

c = 3e5

def SN(flux,eflux,sn):
	A = np.divide(flux,eflux)
	A[A<sn] = np.nan
	return np.divide(A,A)


def rotcur(galaxy, vel_map,evel_map,SN,VSYS,PA,INC,X0,Y0,n_it =5, pixel_scale=1,vary_PA=True,vary_INC=True,vary_XC=True,vary_YC=True,vary_VSYS=True,vary_PHI=True,delta=1, rstart=2, rfinal=50, ring_space=2, frac_pixel=2/3., r_back=0, r_bar_max=20,vmode="circular", save_plots = 1, interpolation = 1, errors = 1,survey = "", config = "", e_ISM = 5 ):
	"""
	vary = [Vrot,Vrad,PA,INC,XC,YC,VSYS,theta,Vtan]
	"""



	vel_map = fits.getdata(vel_map)


	[ny,nx] = vel_map.shape
	ext = [nx/2., -nx/2,-ny/2.,ny/2.]

	
	if evel_map == "":
		
		evel_map = np.ones((ny,nx))
		SN = 1e5
	else:
		evel_map = fits.getdata(evel_map)
		evel_map[evel_map>SN]=np.nan

	evel_map[evel_map>SN]=np.nan

	mask_vel=np.divide(evel_map,evel_map)
	vel_ha = vel_map*mask_vel
	vel_ha[vel_ha==0]=np.nan

	if VSYS == "":
		XC,YC,VSYS,e_vsys = KC(vel_map,X0,Y0,pixel_scale)



	PHI_B = 45
	# List of guess values
	if config == "":
		guess = [50,0,PA,INC,X0,Y0,VSYS,PHI_B]
	else:
		guess = [50,0]

		for res in config_file(config):
				param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), eval(res["val"]), eval(res["min"]), eval(res["max"]) 
				guess.append(val)


	#print(guess, frac_pixel)
	#quit()

	# Starting messenge
	from start_messenge import start
	start(galaxy,guess,vmode,config)



	# Maybe apply a filter to remove isolates pixels in the velocity map?
	vel_ha = filter_isolated_cells(vel_ha,struct=np.ones((3,3)))

	def run_code():


		vary = [True,True,vary_PA,vary_INC,vary_XC,vary_YC,vary_VSYS,vary_PHI,True]
		sigma = []
		if vmode == "bisymmetric":
			PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL,PA_BAR_MAJOR,PA_BAR_MINOR,CHISQ_r,res_MCMC = init_models(vel_ha,evel_map,guess,vary, n_it, sigma, delta, rstart, rfinal, ring_space, frac_pixel,r_back, r_bar_max,pixel_scale,vmode,errors, config, e_ISM)
		else:

			PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL,CHISQ_r,res_MCMC = init_models(vel_ha,evel_map,guess,vary, n_it, sigma, delta, rstart, rfinal, ring_space, frac_pixel,r_back, r_bar_max,pixel_scale,vmode,errors, config, e_ISM)
			PA_BAR_MAJOR,PA_BAR_MINOR = 0,0

		return PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL,PA_BAR_MAJOR,PA_BAR_MINOR,CHISQ_r,res_MCMC




	PA,INC,XC,YC,VSYS,THETA,R,Vrot,Vrad,Vtan,MODEL,PA_BAR_MAJOR,PA_BAR_MINOR,CHISQ_r,res_MCMC = run_code()

	#
	# MCMC to estimate errors on the derived parameters
	#

	#MCMC_results
	
	Vrot_MCMC,e_Vrot,Vrad_MCMC,e_Vrad,PA_MCMC,e_PA,INC_MCMC,e_INC,X0_MCMC,XC_e,YC_MCMC,YC_e,Vsys_MCMC,e_Vsys,theta,e_theta,Vtan_MCMCM,e_Vtan = res_MCMC	

 



	#
	## Write output into a table
	#
	from write_table import write
	if vmode == "circular" or vmode == "radial":
		if survey != "":
			kin_params_table = "ana_kin_%s_model.%s.csv"%(vmode,survey)
		else:
			kin_params_table = "ana_kin_%s_model.%s.csv"%(vmode,galaxy)

		# write header of table
		if path.exists(kin_params_table) == True:
			pass
		else:
			hdr = ["object", "XC", "e_XC_l", "e_XC_u", "YC", "e_YC_l", "e_YC_u", "PA", "e_PA_l", "e_PA_u", "INC", "e_INC_l", "e_INC_u", "VSYS", "e_VSYS_l", "e_VSYS_u", "CHISQ_r" ]
			write(hdr,kin_params_table,column = False)

		if errors == 0:

			kin_params = [galaxy,XC,0,0,YC,0,0,PA,0,0,INC,0,0,VSYS,0,0,CHISQ_r]



		else:

			kin_params = [galaxy,XC,XC_e[0],XC_e[1],YC,YC_e[0],YC_e[1],PA,e_PA[0],e_PA[1],INC,e_INC[0],e_INC[1],VSYS,e_Vsys[0],e_Vsys[1],CHISQ_r]

		write(kin_params,kin_params_table,column = False)

	if vmode == "bisymmetric":
		if survey != "":
			kin_params_table = "ana_kin_%s_model.%s.csv"%(vmode,survey)
		else:
			kin_params_table = "ana_kin_%s_model.%s.csv"%(vmode,galaxy)


		# write header of table
		if path.exists(kin_params_table) == True:
			pass
		else:
			hdr = ["object", "XC", "e_XC_l", "e_XC_u", "YC", "e_YC_l", "e_YC_u", "PA", "e_PA_l", "e_PA_u", "INC", "e_INC_l", "e_INC_u", "VSYS", "e_VSYS_l", "e_VSYS_u", "PHI_BAR", "e_PHI_BAR_l","e_PHI_BAR_u","PA_BAR_MAJOR_sky","PA_BAR_MINOR_sky","CHISQ_r" ]
			write(hdr,kin_params_table,column = False)

		if errors == 0:

			kin_params = [galaxy,XC,0,0,YC,0,0,PA,0,0,INC,0,0,VSYS,0,0,THETA,0,0,PA_BAR_MAJOR,PA_BAR_MINOR,CHISQ_r]



		else:

			kin_params = [galaxy,XC,XC_e[0],XC_e[1],YC,YC_e[0],YC_e[1],PA,e_PA[0],e_PA[1],INC,e_INC[0],e_INC[1],VSYS,e_Vsys[0],e_Vsys[1],THETA,e_theta[0],e_theta[1],PA_BAR_MAJOR,PA_BAR_MINOR,CHISQ_r]

		write(kin_params,kin_params_table,column = False)


	save_file = save_plots


	if interpolation == 1:
		MODEL2D = MODEL[0]
	else:
		MODEL2D = MODEL[1]		

	from plot_models import plot_kin_models
	plot_kin_models(galaxy,vmode,vel_ha,R,Vrot,e_Vrot,Vrad,e_Vrad,Vtan,e_Vtan, VSYS, MODEL2D, ext,plot = 0, save = save_file)



	from kin_model_fits import save_model
	if vmode == "circular": VMODELS = [Vrot]
	if vmode == "radial": VMODELS = [Vrot,Vrad]
	if vmode == "bisymmetric": VMODELS = [Vrot,Vrad,Vtan]

	s = save_model(galaxy,vmode,R,VMODELS,PA,INC,XC,YC,VSYS,save=save_file)


	from vlos_final_model import save_vlos_model
	save_vlos_model(galaxy,vmode,vel_map,MODEL2D,PA,INC,XC,YC,VSYS,save = save_file)


	from fit_rotation_curve_fits import fit_rotcur
	Vmax,Rturn = fit_rotcur(galaxy,vmode,e_Vrot,survey)


	print("Done!")
	print("############################")




 
