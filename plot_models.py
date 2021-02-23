import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
from axis import AXIS 
from CBAR import colorbar as cb
from colormaps_CLC import vel_map
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from error_bars import error_bar

cmap = vel_map()

def plot_kin_models(galaxy,vmode,vel_ha,R,Vrot,eVrot,Vrad,eVrad,Vtan,eVtan,VSYS, MODEL, ext,plot = 0, save = 1):
	e1,e2,e3 = len(eVrot[0]),len(eVrad[0]),len(eVtan[0])

	mask_MODEL = np.divide(MODEL,MODEL)
	MODEL = MODEL*np.divide(vel_ha, vel_ha)
	fig=plt.figure(figsize=(6,2))
	gs2 = GridSpec(1, 3)
	gs2.update(left=0.045, right=0.62,top=0.815,hspace=0.01,bottom=0.135,wspace=0.)


	ax=plt.subplot(gs2[0,0])
	ax1=plt.subplot(gs2[0,1])
	ax2=plt.subplot(gs2[0,2])





	#cmap = plt.cm.get_cmap('bwr')
	#cmap = c_map.reversed()

	vmin = abs(np.nanmin(MODEL))
	vmax = abs(np.nanmax(MODEL))
	max_vel = np.nanmax([vmin,vmax])
	
	vmin = -(max_vel//50 + 1)*50
	vmax = (max_vel//50 + 1)*50


	im0 = ax.imshow(vel_ha - VSYS,cmap = cmap, origin = "lower",vmin = vmin,vmax = vmax, aspect = "auto", extent = ext, interpolation = "nearest")
	im2 = ax1.imshow(MODEL,cmap = cmap, origin = "lower", aspect = "auto", vmin = vmin,vmax = vmax, extent = ext, interpolation = "nearest")

	residual = (vel_ha*mask_MODEL - VSYS)- MODEL
	im2 = ax2.imshow(residual,cmap = cmap, origin = "lower", aspect = "auto",vmin = -50,vmax = 50, extent = ext, interpolation = "nearest")


	AXIS(ax,tickscolor = "k")
	AXIS(ax1,tickscolor = "k",remove_yticks= True)
	AXIS(ax2,tickscolor = "k",remove_yticks= True)



	ax.set_ylabel('$\mathrm{ \Delta Dec~(pix)}$',fontsize=8,labelpad=0)
	ax.set_xlabel('$\mathrm{ \Delta RA~(pix)}$',fontsize=8,labelpad=0)
	ax1.set_xlabel('$\mathrm{ \Delta RA~(pix)}$',fontsize=8,labelpad=0)
	ax2.set_xlabel('$\mathrm{ \Delta RA~(pix)}$',fontsize=8,labelpad=0)

	ax.text(0.05,1.01, "$\mathrm{vlos}$", fontsize = 7, transform = ax.transAxes)
	ax1.text(0.05,1.01,"$\mathrm{model}$", fontsize = 7, transform = ax1.transAxes)
	ax2.text(0.05,1.01,"$\mathrm{residual}$",fontsize = 7, transform = ax2.transAxes)


	#ax.set_facecolor('#e8ebf2')
	#ax1.set_facecolor('#e8ebf2')
	#ax2.set_facecolor('#e8ebf2')


	gs2 = GridSpec(1, 1)
	gs2.update(left=0.68, right=0.995,top=0.815,bottom=0.135)
	ax3=plt.subplot(gs2[0,0])

	from plot_real_mode import plot_true
	#plot_true(vmode,ax3)

	if vmode == "circular":
		ax3.plot(R,Vrot, color = "k",linestyle='-', alpha = 0.6, linewidth=0.8, label = "$\mathrm{V_{t}}$")

		if e1 > 0:
			ax3.errorbar(R,Vrot, yerr=[eVrot[0],eVrot[1]], fmt='o', color = "k",markersize = 1, elinewidth = 0.6,capsize  = 1.5)
			ax3.fill_between(R, Vrot-eVrot[0], Vrot+eVrot[1], color = "darkgray", alpha = 0.4)
			#error_bar(ax3,R,Vrot, eVrot[0],eVrot[1],color = "y")

	if vmode == "radial":

		ax3.plot(R,Vrot, color = "k",linestyle='-', alpha = 0.6, linewidth=0.8, label = "$\mathrm{V_{t}}$")

		if e1 > 0:
			ax3.errorbar(R,Vrot, yerr=[eVrot[0],eVrot[1]], fmt='o', color = "k",markersize = 1.5, elinewidth = 0.6,capsize  = 1.5)
			ax3.fill_between(R, Vrot-eVrot[0], Vrot+eVrot[1], color = "darkgray", alpha = 0.4)
			#error_bar(ax3,R,Vrot, eVrot[0],eVrot[1],color = "darkgray")


		ax3.plot(R,Vrad, color = "orange",linestyle='-', alpha = 0.6, linewidth=0.8, label = "$\mathrm{V_{r}}$")
		if e2 > 0:
			ax3.errorbar(R,Vrad, yerr=[eVrad[0],eVrad[1]], fmt='o', color = "k",markersize = 1.5, elinewidth = 0.6,capsize  = 1.5)
			ax3.fill_between(R, Vrad-eVrad[0], Vrad+eVrad[1], color = "orange", alpha = 0.4)
			#error_bar(ax3,R,Vrad, eVrad[0],eVrad[1],color = "darkorange")




	if vmode == "bisymmetric":
		ax3.plot(R,Vrot, color = "k",linestyle='-', alpha = 0.6, linewidth=0.8, label = "$\mathrm{V_{t}}$")

		if e1 > 0:
			ax3.errorbar(R,Vrot, yerr=[eVrot[0],eVrot[1]], fmt='o', color = "k",markersize = 1.5, elinewidth = 0.6,capsize  = 1.5)
			ax3.fill_between(R, Vrot-eVrot[0], Vrot+eVrot[1], color = "darkgray", alpha = 0.4)


		ax3.plot(R,Vrad, color = "orange",linestyle='-', alpha = 0.6, linewidth=0.8, label = "$\mathrm{V_{2,r}}$")

		if e2 > 0:
			ax3.errorbar(R,Vrad, yerr=[eVrad[0],eVrad[1]], fmt='o', color = "orange",markersize = 1.5, elinewidth = 0.6,capsize  = 1.5)
			ax3.fill_between(R, Vrad-eVrad[0], Vrad+eVrad[1], color = "orange", alpha = 0.4)




		ax3.plot(R,Vtan, color = "skyblue",linestyle='-', alpha = 0.6, linewidth=0.8, label = "$\mathrm{V_{2,t}}$")

		if e3 > 0:
			ax3.errorbar(R,Vtan, yerr=[eVtan[0],eVtan[1]], fmt='o', color = "dodgerblue",markersize = 1.5, elinewidth = 0.6,capsize  = 1.5)
			ax3.fill_between(R, Vtan-eVtan[0], Vtan+eVtan[1], color = "dodgerblue", alpha = 0.4)





	ax3.legend(loc = "center", fontsize = 6.5, bbox_to_anchor = (0, 1, 1, 0.1), ncol = 3, frameon = False)

	vels = [Vrot, Vrad, Vtan]
	max_vel,min_vel = int(np.nanmax(vels)),int(np.nanmin(vels)) 

	ax3.set_ylim(min_vel-50,max_vel+50)
	ax3.plot([0,np.nanmax(R)],[0,0],color = "k",linestyle='-', alpha = 0.6,linewidth = 0.5)
	ax3.set_xlabel('$\mathrm{r~(arcsec)}$',fontsize=8,labelpad = 0)
	ax3.set_ylabel('$\mathrm{V_{rot}~(km~s^{-1})}$',fontsize=8,labelpad = 0)
	ax3.set_facecolor('#e8ebf2')

	AXIS(ax3,tickscolor = "k")

	# Make a plot with major ticks that are multiples of 10 and minor ticks that
	# are multiples of 5.  Label major ticks with '%d' formatting but don't label
	# minor ticks.
	if np.nanmax(R) // 10 > 1:
		ax3.xaxis.set_major_locator(MultipleLocator(10))
		ax3.xaxis.set_major_formatter(FormatStrFormatter('%d'))
		# For the minor ticks, use no labels; default NullFormatter.
		ax3.xaxis.set_minor_locator(MultipleLocator(2))
	else:
		ax3.xaxis.set_major_locator(MultipleLocator(2))
		ax3.xaxis.set_major_formatter(FormatStrFormatter('%d'))
		# For the minor ticks, use no labels; default NullFormatter.
		ax3.xaxis.set_minor_locator(MultipleLocator(1))


	if np.nanmax(Vrot) // 100 > 1:
		ax3.yaxis.set_major_locator(MultipleLocator(100))
		ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))
		# For the minor ticks, use no labels; default NullFormatter.
		ax3.yaxis.set_minor_locator(MultipleLocator(20))
	else:
		ax3.yaxis.set_major_locator(MultipleLocator(50))
		ax3.yaxis.set_major_formatter(FormatStrFormatter('%d'))
		# For the minor ticks, use no labels; default NullFormatter.
		ax3.yaxis.set_minor_locator(MultipleLocator(25))



	cb(im0,ax,orientation = "horizontal", colormap = cmap, bbox= (0.5,1.135,1,1),width = "100%", height = "5%",label_pad = -17, label = "$\mathrm{(km~s^{-1})}$",font_size=6)
	cb(im2,ax2,orientation = "horizontal", colormap = cmap, bbox= (0,1.135,1,1),width = "100%", height = "5%",label_pad = -17, label = "$\mathrm{(km~s^{-1})}$",font_size=6)


	if save == 1 and plot == 1:
		plt.savefig("./plots/kin_%s_model_%s.png"%(vmode,galaxy),dpi = 300)
		plt.show()
		plt.clf()
	else:

		if plot == 1:
			plt.show()
			plt.clf()
		if save == 1:
			plt.savefig("./plots/kin_%s_model_%s.png"%(vmode,galaxy),dpi = 300)
			plt.clf()




