
Installation

(1) Make a directory to put the executables in:
mkdir ~/XS-code

(2) copy all files in this new folder

(3) Make the file executable:
chmod +x main_xs.py


(4) Edit  .bashrc:
gedit ~/.bashrc

Add the folowing line:
export PATH="$PATH:$HOME/XS-code"

Refresh the bash file:
source ~/.bashrc

Try it. Go to any folder and type  main_xs




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
residual map


* plots
plot of the best 2D kinematic models
plot of the asymptotic velocity estimated with Vt 



************************************************
************    ADDITIONAL *********************
************************************************
The following directories need to be created in the working directory
./plots
./fits
./vmax_rturn







#
# RUNNING THE CODE:
#
#
# EXAMPLE
#

main_xs.py test_galaxy test.fits "" 10 "" 80 50 61 25 5 1 1 1 1 1 1 1 2.5 "" 3 60 "" 1/3. 0 120 circular 1 1 0 "" "" ""

