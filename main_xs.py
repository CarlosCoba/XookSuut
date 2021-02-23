#! /usr/bin/python3
import sys
nargs=len(sys.argv)
from initialize_XS import rotcur

"""
#################################################################################
# 				The XookSuut-code.				#
#				version 23.02.17				#
# 				C. Lopez-Coba					#
#################################################################################


INPUT:
object = 	object name.
vel_map = 	2D velocity map in km/s.
evel_map = 	2D error map in km/s. if no error map then set this to ""
SN = 		Maximum value in the error map allowed for computing the rotation curve.
		Better results are obtained when SN has low values 15-20 km/s.
		if evel_map = 0 then the whole velocity map will be used for estimating the rotation curve.
VSYS	=	Guess value for Vsys. If VSYS = "", then it will be used the central value given by VSYS =  vel_map[Y0][X0].
PA	= 	Guess value for the Position Angle of the major axis.		
INC	=	Guess value for the Inclination.
X0,Y0	=	Guess values for the kinematic center. 
N_it 	=	Number of iteratios to find the best parameters. Set to 5.
pixel_scale	Size of the pixel in arcsecs.
vary_PA =	Boolean.
vary_INC =	Boolean.
vary_XC =	Boolean.
vary_YC =	Boolean.
vary_VSYS =	Boolean.
vary_PHI_bar =	Boolean.
FWHM	=	PSF FWHM. If FWHM is given then Delta is set to FWHM/2. and ring_space = FWHM. This ensure no free space between rings.
Delta =		width_ring =	2*Delta [arcsec]. If FWHM is given, then this value is not taken into account.
ring_space =	the spacing between rings [arcsec]. If FWHM is given, then this value is not taken into account.
rstart	=	the radii of the annulis were the fitting is going to start [arcsec].
rfinal	=	Maximum radius where the kinematic model will be computed [arcsec].
frac_pixel =	The fraction of pixels within each annuli needed to calculate the model in each ring.
r_back	=	If you want to explore backwards rings. [arcsec].
r_bar_max= 	Maximum lengh of the bar like perturbation [arcsec].
model	= 	The different kinematic models. It can be pure circular motion [circular], radial flows [radial], or barlike flows [bisymmetric].
errors	=	Boolean. If you want to compute the errors via MCMC for the derived parameters.
survey =	String. If the object belongs to a specific galaxy survey.
config =	configure file to pass initial guess and constrains for the constant params (x0, y0, Vsys, pa, inc, phi_bar). This file is optional.
		If this file is passed to main_xs, then it ommits the previos guess values (VSYS, PA, INC, X0, Y0) as well as the VARY_ entrance

save_plots =	boolean
interp = 	boolean. This performs a linear interpolation in the best 2D model.
e_ISM =		Error in the emission line centroid.


*** Values in brackets are optional

A configuration file has the follow entrance

#
#
# XS config file
# Example of configuration file
# param col. are the constant parameters to fit. Do not touch this column.
# val col. is the initial values for the considered parameters
# fit col. is a boolen. If set to 1 means the parameter is fitted, other wise set 0.
# min col. is the minimum value for the considered parameter. If fit is seto to 0, the min value is not taken into account. 
# max col. is the maximum value for the considered parameter. If fit is seto to 0, the max value is not taken into account. 
#
#
#
param	val	fit	min	max
pa	35	1	0	360
inc	35	1	0	90
x0	25	1	0	50
y0	25	1	0	50
Vsys	11168.266579019075	0	0	3e6
phi_b	45	1	0	180





					OUTPUT:


-Tables:

*Table containing the best fitted values for the constant parameters (x0, y0, Vsys, pa, inc, phi_bar)
ana_kin_model.csv

*Table containing different estimations of the maximum circular velocity.
vmax_rc.model.csv 


* fits files
2D array of the LoV of the best kinematic model
1D array of the different kinematic models as function of the deprojected radius


* plots
plot of the best 2D kinematic models
plot of the asymptotic velocity estimated with Vt 


ADDITIONAL
The following directories need to be created in the same path
./plots
./fits
./vmax_rturn
"""


if __name__ == "__main__":


	if (nargs == 33 ):
		galaxy = sys.argv[1]
		vel_map = sys.argv[2]
		evel_map = sys.argv[3]
		SN = float(sys.argv[4])
		VSYS = sys.argv[5]
		PA = float(sys.argv[6])
		INC = float(sys.argv[7])
		X0 = float(sys.argv[8])
		Y0 = float(sys.argv[9])
		n_it = sys.argv[10]
		pixel_scale = float(sys.argv[11])
		vary_PA,vary_INC,vary_XC,vary_YC,vary_VSYS, vary_PHI = bool(float(sys.argv[12])), bool(float(sys.argv[13])), bool(float(sys.argv[14])), bool(float(sys.argv[15])), 			bool(float(sys.argv[16])), bool(float(sys.argv[17]))
		FWHM_arc = float(sys.argv[18])
		delta = sys.argv[19]


		rstart = float(sys.argv[20])
		rfinal = float(sys.argv[21])
		ring_space = sys.argv[22]


		frac_pixel = eval(sys.argv[23])
		r_back = float(sys.argv[24])
		r_bar_max = float(sys.argv[25])
		vmode = sys.argv[26]

		save_plots = bool(float(sys.argv[27]))
		interpolation = bool(float(sys.argv[28]))
		errors = bool(float(sys.argv[29]))
		survey = sys.argv[30]
		config = sys.argv[31]
		e_ISM = sys.argv[32]

		if VSYS != "": float(sys.argv[5])
		if n_it == "": n_it = 5
		if n_it != "": n_it = int(float(n_it))
		if FWHM_arc != "":
			delta = FWHM_arc/2. 
			ring_space = FWHM_arc


		if pixel_scale == "": pixe_scale = 1
		if rstart == "": rstart = 2.5
		if rfinal == "": rfinal = 40
		if frac_pixel == "": frac_pixel = 2/3.
		if r_back == "": r_back = 0
		if r_bar_max == "": r_bar_max = 1e3
		if errors == "": errors = bool(0)
		if survey == "": survey = ""
		if config == "": config = ""
		if save_plots == "": save_plots = 1
		if interpolation == "": interpolation = 1
		if e_ISM == "": e_ISM = 5
		if e_ISM != "": e_ISM = float(e_ISM)


		if delta != "": delta = float(delta)
		if ring_space != "": ring_space = float(ring_space)


		rotcur(galaxy, vel_map, evel_map, SN, VSYS, PA, INC, X0, Y0, n_it, pixel_scale, vary_PA, vary_INC, vary_XC, vary_YC, vary_VSYS, vary_PHI, delta, rstart, rfinal, ring_space, frac_pixel, 		r_back, r_bar_max,vmode, save_plots, interpolation, errors, survey,config,e_ISM)


	else:
		print ("USE: XS-code.py object vel_map [evel_map] SN [VSYS] PA INC X0 Y0 [N_it=5] [pixel_scale=1] vary_PA[0,1] vary_INC[0,1] vary_XC[0,1] vary_YC[0,1] vary_VSYS[0,1] vary_THETA_bar[0,1] [FWHM_arc] Delta [Rstart=2.5] [Rfinal=40] spacing [frac_pixel=2/3.] [R_back=0] [R_bar_max=1e3] [model=cicular,radial,bisymmetric] [save_plots=1] [interp=1] [errors=0] [survey] [config_file] [e_ISM=5]")

		exit()


"""


#
# RUNNING THE CODE:
#
#
# EXAMPLE
#

main_xs.py test_galaxy test.fits "" 10 "" 80 50 61 25 5 1 1 1 1 1 1 1 2.5 "" 3 60 "" 1/3. 0 120 circular 1 1 0 "" "" ""


"""





