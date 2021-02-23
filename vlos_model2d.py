import numpy as np



def vlos(x,y,Vrot,Vr2,pa,inc,x0,y0,Vsys,Vt2,phi_b_sky,vmode):
	pa,inc,phi_b_sky=(pa)*np.pi/180,inc*np.pi/180,phi_b_sky*np.pi/180

	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))
	R= np.sqrt(X**2+(Y/np.cos(inc))**2)

	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)
	m = 2


	theta = np.arctan(sin_tetha/cos_tetha)
	#phi_b = np.arctan(np.tan(phi_b_sky-pa)/np.cos(inc))
	phi_b = phi_b_sky
	theta_b = theta - phi_b


	if vmode == "circular":
		vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha)

	if vmode == "radial":
		vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha + Vr2*sin_tetha)

	if vmode == "bisymmetric":	
		vlos = Vsys+np.sin(inc)*((Vrot*cos_tetha)-Vt2*np.cos(m*theta_b)*cos_tetha - Vr2*np.sin(m*theta_b)*sin_tetha)


	return np.ravel(vlos)


