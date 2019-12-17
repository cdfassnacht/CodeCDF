import numpy as np
import matplotlib.pylab as plt


#  function that calculates the Power spectrum for a field, where the pixels correspond to kabs, and the kd determines the bins in k - space
def get_power(field,kabs, kd):
	kabsflat = kabs.flatten()
	fft_field_abs = np.absolute( np.fft.fft2(field) )**2.0
	
	aux0 = ( fft_field_abs ).flatten()
	aux1, k = np.histogram(kabsflat,bins = kd)
	aux2, _  = np.histogram(kabsflat,weights = aux0,bins=kd)
	aux3, _  = np.histogram(kabsflat,weights = aux0**2,bins=kd)
	aux1[np.where(aux1==0)] += 1
	return  aux2/(aux1) / field.shape[0]/ field.shape[1] ,  ( k[1:] + k[:-1] )/2.0 ,  np.sqrt( aux3/(aux1) -  ( aux2/(aux1) )**2 ) / field.shape[0]/ field.shape[1]
	
#  get the xy grid, kxky grid for a given box
def get_grids(xL,yL,xN,yN):

	x = np.linspace(-xL/2.0,xL/2.0,xN)
	y = np.linspace(-xL/2.0,xL/2.0,xN)

	xk = np.fft.fftfreq(xN)* np.pi*2.0 /xL
	yk = np.fft.fftfreq(yN)* np.pi*2.0 /yL

	xy =  np.array([np.ones( (yN))[:,np.newaxis]*x[np.newaxis,:],np.ones((xN))[np.newaxis,:]*y[:,np.newaxis]])
	
	kxky =  np.array([np.ones( (yN))[:,np.newaxis]*xk[np.newaxis,:],np.ones((xN))[np.newaxis,:]*yk[:,np.newaxis]])
	

	return xy , kxky
	


# generate fields	
N = 1000
L = 1.0

xy, kxky = get_grids(L,L,N,N)	
	
kabs =  np.sqrt(np.sum(kxky**2,axis=0))

K=100 
kd = np.linspace(np.min(kabs),np.max(kabs),K)
white_noise = np.random.normal(0.0,1.0,(N,N))
power_field_analytical = lambda k :  ( np.power(k+0.1,-4)+0.5 ) /10.0
white_noise_fft = np.fft.fft2(white_noise)  * np.sqrt(  power_field_analytical(kabs) )

field =  np.fft.ifft2(white_noise_fft).real



# get the power spectrum, mean value of the k-bins, the standard deviation of the power spectrum
P,k,P_std = get_power(white_noise,kabs, kd)
P_field,k_field,P_std_field = get_power(field,kabs, kd)




# plotting

plt.subplot(121)
plt.plot(k_field,power_field_analytical(k_field),c="r")
plt.plot(k,P,"b")
#plt.plot(k,P-P_std,ls="--",c="b")
#plt.plot(k,P+P_std,ls="--",c="b")
plt.plot(k_field,P_field,c="g",ls="--")
#plt.plot(k_field,P_field-P_std_field,ls="--",c="g")
#plt.plot(k_field,P_field+P_std_field,ls="--",c="g")
plt.yscale("log")
plt.xscale("log")
plt.subplot(122)
plt.imshow(field)

plt.show()
