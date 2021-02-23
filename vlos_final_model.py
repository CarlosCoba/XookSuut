import numpy as np
from astropy.io import fits



def save_vlos_model(galaxy,vmode,vel_map,VLOS,PA,INC,XC,YC,VSYS,save = 1):

	if save == 1:

		# The best 2D model
		hdu = fits.ImageHDU()
		hdu.data = VLOS

	
		hdu.header['NAME0'] = '2D LoS velocity model'
    
		hdu.header['PA'] = PA
		hdu.header['INC'] = INC
		hdu.header['VSYS'] = VSYS
		hdu.header['XC'] = XC
		hdu.header['YC'] = YC
		
		hdu.writeto("./fits/%s.%s.vlos_model.fits"%(galaxy,vmode),overwrite=True)


		
		# Now the residual map

		hdu = fits.ImageHDU()
		hdu.data = vel_map-(VLOS+VSYS)

	
		hdu.header['NAME0'] = 'residual map, data-model'
    
		hdu.header['PA'] = PA
		hdu.header['INC'] = INC
		hdu.header['VSYS'] = VSYS
		hdu.header['XC'] = XC
		hdu.header['YC'] = YC
		
		hdu.writeto("./fits/%s.%s.residual.fits"%(galaxy,vmode),overwrite=True)



	else: pass

