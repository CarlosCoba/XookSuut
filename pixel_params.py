import numpy as np
import numpy as np
import matplotlib.pylab as plt
import sys
import lmfit
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
from matplotlib.gridspec import GridSpec


 
def Rings(xy_mesh,pa,inc,x0,y0,pixel_scale):
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	(x,y) = xy_mesh
	#X = -(x-x0)*np.sin(PA)+(y-y0)*np.cos(PA)
	#Y = ((x-x0)*np.cos(PA)+(y-y0)*np.sin(PA))/np.cos(inc)

	X = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))
	Y = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))



	R= np.sqrt(X**2+(Y/np.cos(inc))**2)
	return R*pixel_scale

def Rings0(xy_mesh,pa,inc,x0,y0):
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	(x,y) = xy_mesh
	X = -(x-x0)*np.sin(PA)+(y-y0)*np.cos(PA)
	Y = ((x-x0)*np.cos(PA)+(y-y0)*np.sin(PA))/np.cos(inc)
	R= np.sqrt(X**2+Y**2)
	return R



def Rings_r_1(R,vel,ring,delta=0):
	[ny,nx]= vel.shape
	M=np.ones((ny,nx))

	if ring == 0:
		mask = (R>=ring) & (R <= ring+delta)
		#mask = R==0
	else:
		#ring = ring +1
		ring = ring +delta
		mask = (R>=ring) & (R <= ring+delta)		
	
	s = M*mask
	return s



def Rings_r_2(R,vel,ring,delta=0,pixel_scale=1):
	R = R*pixel_scale
	[ny,nx]= vel.shape
	M=np.ones((ny,nx))

	if ring == 0:
		mask = (R>=ring) & (R <= ring+delta)
		#mask = R==0
	else:
		#ring = ring +1
		ring = ring*delta +delta
		mask = (R>=ring) & (R <= ring+delta)		
	
	s = M*mask
	#plt.imshow(s)
	#plt.show()
	return s


def Rings_r_3(R,vel,ring,delta=0):
	[ny,nx]= vel.shape
	M=np.ones((ny,nx))

	if ring == 0:
		mask = (R>=ring) & (R <= ring+delta)
	else:
		mask = (R>=ring) & (R <= ring+delta)		
	
	s = M*mask
	return s


def Rings_r_4(R,vel,ring,delta=0):
	[ny,nx]= vel.shape
	M=np.ones((ny,nx))

	if ring == 0:
		mask = (R>=ring) & (R <= ring+delta)
	else:
		mask = (R>=ring-0.5*delta) & (R <= ring+0.5*delta)		
	
	s = M*mask
	#plt.imshow(s)
	#plt.show()
	return s


def Ring_arcsec(R,vel,ring,delta=0):
	R = R
	if delta == 0:
		R = R.astype(int)


	[ny,nx]= vel.shape
	M=np.ones((ny,nx))

	if delta == 0:
		R = R.astype(int)
		mask = R== ring

	else:

		mask = (R>=ring-delta) & (R < ring+delta)



	s = M*mask

	return s



def Vlos_BISYM(xy_mesh,Vrot,Vr2,pa,inc,x0,y0,Vsys):
	(x,y) = xy_mesh
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	cos_tetha = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))/R
	sin_tetha = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))/(np.cos(inc)*R)
	#vlos =  Vrot*cos_tetha*np.sin(inc)+Vsys
	#m = 2
	#vlos = Vsys+np.sin(inc)*((Vrot*cos_tetha)-Vt2*np.cos(m*theta_b*np.pi/180)*cos_tetha - Vr2*np.sin(m*theta_b*np.pi/180)*sin_tetha)
	vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha + Vr2*sin_tetha)
	return np.ravel(vlos)


def Vlos_ROT(xy_mesh,Vrot,pa,inc,x0,y0,Vsys):#,Vt2=0,Vr2=0,theta_b=0):
	(x,y) = xy_mesh
	R  = Rings0(xy_mesh,pa,inc,x0,y0)
	#print "R=",R
	PA,inc=(pa-90*0)*np.pi/180,inc*np.pi/180
	cos_tetha = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))/R
	sin_tetha = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))/(np.cos(inc)*R)
	vlos =  Vrot*cos_tetha*np.sin(inc)+Vsys
	return np.ravel(vlos)



def Vrot_MODEL(xy_mesh,vlos,pa,inc,x0,y0,Vsys):#,Vt2=0,Vr2=0,theta_b=0):
	(x,y) = xy_mesh
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	PA,inc=(pa-90*0)*np.pi/180,inc*np.pi/180
	cos_tetha = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))/R
	sin_tetha = (- (x-x0)*np.cos(PA) - (y-y0)*np.sin(PA))/(np.cos(inc)*R)
	vrot = (vlos - Vsys)/cos_tetha*np.sin(inc)
	return vrot


def Weight(xy_mesh,Vrot,pa,inc,x0,y0,Vsys):
	(x,y) = xy_mesh
	PA,inc=(pa)*np.pi/180,inc*np.pi/180
	R  = Rings0(xy_mesh,pa,inc,x0,y0)
	cos_tetha = (- (x-x0)*np.sin(PA) + (y-y0)*np.cos(PA))/R
	abs_cos = abs(cos_tetha) 
	return np.ravel(abs_cos)







def pixels(vel,evel,guess,ringpos, delta=1,ring = "pixel",pixel_scale = 1):

	vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0 = guess
	[ny,nx] = vel.shape

	x = np.linspace(0, nx, nx)
	y = np.linspace(0, ny, ny)

	x = np.arange(0, nx, 1)
	y = np.arange(0, ny, 1)

	XY_mesh = np.meshgrid(x,y,sparse=True)
	indices = np.indices((ny,nx))

	R_galaxy_plane = Rings(XY_mesh,guess[2],guess[3],guess[4],guess[5],pixel_scale)
	R_galaxy_plane = R_galaxy_plane*np.divide(vel,vel)


	pix_y=  indices[0]
	pix_x=  indices[1]

	R_galaxy = Rings(XY_mesh,guess[2],guess[3],guess[4],guess[5],pixel_scale)
	R_int= R_galaxy.astype(int)

	i = ringpos


	if ring == "arcsec" : 
		R= Rings_r_2(R_int,vel,i,delta,pixel_scale)
	if ring == "pixel" : 
		R= Rings_r_4(R_galaxy,vel,i,delta)
	if ring == "broad":
		R= Rings_r_3(R_int,vel,i,delta)
	if ring == "ARCSEC":
		#R= Ring_arcsec(R_galaxy,vel,int(i),delta)
		R= Ring_arcsec(R_galaxy,vel,i,delta)
		#R= Ring_arcsec(R_int,vel,int(i),delta,pixel_scale)



	mask = R ==1
	indices = np.indices((ny,nx))


	pix_y=  indices[0]
	pix_x=  indices[1]


	pix_x = pix_x[mask]
	pix_y = pix_y[mask]




	vel_val = [vel[m][n] for m,n in zip(pix_y,pix_x)]
	R_val = [R_int[m][n] for m,n in zip(pix_y,pix_x)]
	#R_val = [R[m][n] for m,n in zip(pix_x,pix_y)]
	e_vel = [evel[m][n] for m,n in zip(pix_y,pix_x)]


	"""
	#Check if there it no abrupt changes in the velocity
	vel_0 = np.nanmedian(vel_val)
	Vel_val=[]

	Pix_x,Pix_y,E_vel,r_val = [],[],[],[]
	for w in range(len(vel_val)):


		if np.isfinite(vel_val[w]) == True:
			diff = abs(vel_val[w]-vel_0)
			if diff < 400:
				#vel_0 = vel_val[w]
				Vel_val.append(vel_val[w])
			
				Pix_x.append(pix_x[w])
				Pix_y.append(pix_y[w])
				E_vel.append(e_vel[w])
				r_val.append(R_val[w])

			else:

				Vel_val.append(np.nan)
				Pix_x.append(pix_x[w])
				Pix_y.append(pix_y[w])
				E_vel.append(e_vel[w])
				r_val.append(R_val[w])

		else:
				#print("aquiiiiiiiiii",e_vel[w])

				Vel_val.append(np.nan)
				Pix_x.append(pix_x[w])
				Pix_y.append(pix_y[w])
				E_vel.append(e_vel[w])
				r_val.append(R_val[w])
		#	Vel_val.append(np.nan)


	vel_val = Vel_val
	pix_x,pix_y,e_vel,R_val= Pix_x,Pix_y,E_vel,r_val
	"""
	
	vel_val = np.array(vel_val)
	R_val = np.array(R_val)
	e_vel = np.array(e_vel)
	pix_y =  np.asarray(pix_y)
	pix_x = np.asarray(pix_x)
	npix_exp = len(pix_x)



	mask = np.isfinite(vel_val) == True
	vel_val = vel_val[mask]
	R_val = R_val[mask]
	e_vel = e_vel[mask]
	pix_y = pix_y[mask]
	pix_x = pix_x[mask]
	#e_vel = 1./e_vel

	len_vel = len(vel_val)


	"""
	plt.imshow(vel, origin = "lower")
	plt.plot(guess[4],guess[5],"gx")
	plt.plot(pix_x,pix_y,"ko")
	plt.show()
	"""
	if npix_exp >0 and len_vel >0 :
		f_pixel = len(vel_val)/(1.0*npix_exp)
	else:
		f_pixel = 0


	return [pix_x,pix_y], vel_val, e_vel, f_pixel







