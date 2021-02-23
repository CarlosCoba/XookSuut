import numpy as np

 
def Rings(xy_mesh,pa,inc,x0,y0,integer = False,pixel_scale=1):
	(x,y) = xy_mesh

	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))

	R= np.sqrt(X**2+(Y/np.cos(inc))**2)
	R = R*pixel_scale

	if integer == True:
		R =R.astype(int)


	return R





def Vlos_BISYM(xy_mesh,Vrot,Vr2,pa,inc,x0,y0,Vsys):
	(x,y) = xy_mesh
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	cos_tetha = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))/R
	sin_tetha = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))/(np.cos(inc)*R)
	vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha + Vr2*sin_tetha)
	return np.ravel(vlos)


def Phi_bar_sky(pa,inc,phi_b_gal):

	pa, inc , phi_b_gal = pa*np.pi/180, inc*np.pi/180 , phi_b_gal*np.pi/180
	phi_sky = pa + np.arctan(np.tan(phi_b_gal)*np.cos(inc))
	return phi_sky*180/np.pi




def weigths_w(xy_mesh,pa,inc,x0,y0,ak,ak_plus_1):

	#where {ak} are the semimajor axes of the K_d ellipses at which
	#the model disk intensities are tabulated

	ai = Rings(xy_mesh,pa,inc,x0,y0)
	delt = (ak_plus_1 - ak)

	w_k_i = (ak_plus_1 - ai) / delta

	w_k_plus_1_i = (ai - ak) / delta


	return w_k_i,w_k_plus_1_i

def cos_sin(xy_mesh,pa,inc,x0,y0,pixel_scale):
	(x,y) = xy_mesh
	pa,inc=(pa)*np.pi/180,inc*np.pi/180
	R  = Rings(xy_mesh,pa,inc,x0,y0,pixel_scale)
	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)
	return np.ravel(cos_tetha),np.ravel(sin_tetha)


def trigonometric_weights(xy_mesh,pa,inc,x0,y0,phi_b_sky,vmode="radial",pixel_scale=1):
	phi_b_sky = phi_b_sky*np.pi/180
	cos,sin = cos_sin(xy_mesh,pa,inc,x0,y0,pixel_scale)


	if vmode == "circular":
		w_rot = np.sin(inc*np.pi/180)*cos
		return w_rot




	if vmode == "radial":
		w_rot = np.sin(inc*np.pi/180)*cos
		w_rad = np.sin(inc*np.pi/180)*sin
		return w_rot,w_rad

	if vmode == "bisymmetric":



		theta = np.arctan(sin/cos)
		#phi_b = np.arctan(np.tan(phi_b_sky-pa*np.pi/180)/np.cos(inc*np.pi/180))
		phi_b = phi_b_sky
		theta_b = theta - phi_b


		w_rot = np.sin(inc*np.pi/180)*cos
		w_rad = np.sin(inc*np.pi/180)*sin*np.sin(2*theta_b)
		w_tan = np.sin(inc*np.pi/180)*cos*np.cos(2*theta_b)
		return w_rot,w_rad,w_tan




def M_radial(xy_mesh,Vrot,Vr2,pa,inc,x0,y0,Vsys,Vt2,theta_b,vel_val,e_vel,vmode = "radial",pixel_scale=1):
	#vmode = "radial"

	if vmode == "circular":

		w_rot = trigonometric_weights(xy_mesh,pa,inc,x0,y0,0,vmode,pixel_scale =pixel_scale)#*weigths_w
		w_rad = 0
		w_sys = 1


		sigma_v = e_vel
		x11,x12 = w_rot**2/sigma_v**2,w_rot*w_rad/sigma_v**2
		x21,x22 = w_rot*w_rad/sigma_v**2,w_rad**2/sigma_v**2


		D = (vel_val)
		y1 = (w_rot/sigma_v**2)*D
		y2 = (w_rad/sigma_v**2)*D

		A = np.asarray([[np.nansum(x11),np.nansum(x12)],[np.nansum(x21),np.nansum(x22)]])
		B= np.asarray([np.nansum(y1),np.nansum(y2)])


		vrot,vrad = np.nansum(y1)/np.nansum(x11), 0

		if np.isfinite(vrot) == False: vrot = 0

		return w_rot,w_rad,vrot,vrad



	if vmode == "radial":
		w_rot = trigonometric_weights(xy_mesh,pa,inc,x0,y0,0,pixel_scale = pixel_scale)[0]#*weigths_w
		w_rad = trigonometric_weights(xy_mesh,pa,inc,x0,y0,0)[1]#*weigths_w



		sigma_v = e_vel
		x11,x12 = w_rot**2/sigma_v**2,w_rot*w_rad/sigma_v**2
		x21,x22 = w_rot*w_rad/sigma_v**2,w_rad**2/sigma_v**2


		D = (vel_val)
		y1 = (w_rot/sigma_v**2)*D
		y2 = (w_rad/sigma_v**2)*D

			

		A = np.asarray([[np.nansum(x11),np.nansum(x12)],[np.nansum(x21),np.nansum(x22)]])
		B= np.asarray([np.nansum(y1),np.nansum(y2)])


		x = np.linalg.solve(A, B)
		vrot,vrad = abs(x[0]),x[1]


		if np.isfinite(vrot) == False: vrot = 0
		if np.isfinite(vrad) == False: vrad = 0

		return w_rot,w_rad,vrot,vrad

	if vmode == "bisymmetric":

		w_rot = trigonometric_weights(xy_mesh,pa,inc,x0,y0,theta_b,vmode)[0]#*weigths_w
		w_rad = trigonometric_weights(xy_mesh,pa,inc,x0,y0,theta_b,vmode)[1]#*weigths_w
		w_tan = trigonometric_weights(xy_mesh,pa,inc,x0,y0,theta_b,vmode)[2]



		sigma_v = e_vel
		x11,x12,x13 = w_rot**2/sigma_v**2,w_rot*w_rad/sigma_v**2,w_tan*w_rot/sigma_v**2
		x21,x22,x23 = w_rot*w_rad/sigma_v**2,w_rad**2/sigma_v**2,w_rad*w_tan/sigma_v**2
		x31,x32,x33 = w_rot*w_tan/sigma_v**2,w_rad*w_tan/sigma_v**2,w_tan**2/sigma_v**2


		D = (vel_val)
		y1 = (w_rot/sigma_v**2)*D
		y2 = (w_rad/sigma_v**2)*D
		y3 = (w_tan/sigma_v**2)*D


		A = np.asarray([[np.nansum(x11),np.nansum(x12),np.nansum(x13)],[np.nansum(x21),np.nansum(x22),np.nansum(x23)],[np.nansum(x31),np.nansum(x32),np.nansum(x33)]])
		B= np.asarray([np.nansum(y1),np.nansum(y2),np.nansum(y3)])



		try:
			x = np.linalg.solve(A, B)
			vrot,vrad,vtan = abs(x[0]),x[1],x[2]
			if np.isfinite(vrot) == False: vrot = 0
			if np.isfinite(vrad) == False: vrad = 0
			if np.isfinite(vtan) == False: vtan = 0

		except(TypeError):
			w_sys,w_rot,w_rad,w_tan,vrot,vrad,vsys,vtan =  0,0,0,0,0,0,0,0

		vsys, w_sys = 0,0
		return w_sys,w_rot,w_rad,w_tan,vrot,vrad,vsys,vtan



