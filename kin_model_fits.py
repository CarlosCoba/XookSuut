import numpy as np
from astropy.io import fits



def save_model(galaxy,vmode,R,MODELS,PA,INC,XC,YC,VSYS,save = 1):
	m = len(MODELS)
	n = len(MODELS[0])

	if vmode == "circular":
			vrot = MODELS[0]
			data = np.zeros((2,n))
			data[0][:] = R
			data[1][:] = vrot

	if vmode == "radial":
			vrot = MODELS[0]
			vrad = MODELS[1]
			data = np.zeros((3,n))
			data[0][:] = R
			data[1][:] = vrot
			data[2][:] = vrad


	if vmode == "bisymmetric":
			vrot = MODELS[0]
			vrad = MODELS[1]
			vtan = MODELS[2]
			data = np.zeros((4,n))
			data[0][:] = R
			data[1][:] = vrot
			data[2][:] = vrad
			data[3][:] = vtan


	if save == 1:
		hdu = fits.ImageHDU()
		hdu.data = data

		if vmode == "circular":
			hdu.header['NAME0'] = 'deprojected distance (arcsec)'
			hdu.header['NAME1'] = 'circular velocity (km/s)     '
		if vmode == "radial":
			hdu.header['NAME0'] = 'deprojected distance (arcsec)'
			hdu.header['NAME1'] = 'circular velocity (km/s)     '
			hdu.header['NAME2'] = 'radial velocity (km/s)     '
		if vmode == "bisymmetric":
			hdu.header['NAME0'] = 'deprojected distance (arcsec)'
			hdu.header['NAME1'] = 'circular velocity (km/s)     '
			hdu.header['NAME2'] = 'radial velocity (km/s)     '
			hdu.header['NAME3'] = 'tangencial velocity (km/s)     '

		hdu.header['PA'] = PA
		hdu.header['INC'] = INC
		hdu.header['VSYS'] = VSYS
		hdu.header['XC'] = XC
		hdu.header['YC'] = YC
		
		hdu.writeto("./fits/%s.%s_model.fits"%(galaxy,vmode),overwrite=True)

