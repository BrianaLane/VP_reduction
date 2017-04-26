from numpy import *
import string
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal 
import pyfits
import subprocess
import numpy.ma as ma
from pylab import *
import os.path as op
import sys
import numpy.ma as ma


#define ID to perform flux calibration on
id = 3
fiberextract = False
create_standSpec = True
create_sensFunc = True

#define wavelength range (in angstroms)
wl_lower = 4700
wl_higher = 7200

#set up path to your cure bin
#extinc_McD.dat needs to be saved in the scripts directory in your CURE folder
curebin_path  = '/Users/Briana/Documents/cure/virusp1/bin/'

#----------#
# Get Data #
#----------#

#define all the files to use and correspoinding error files
if id == 1:
	name = 'Feige34'
	#this should be the file in AB magnitudes
	stdfile= '/Users/Briana/Documents/cure/virusp1/scripts/standard_stars/feige_34_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=123
	#define all the files to use and correspoinding error files
elif id == 2:
	name = 'Feige67'
	#this should be the file in AB magnitudes
	stdfile= '/Users/Briana/Documents/cure/virusp1/scripts/standard_stars/feige_67_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=123
elif id == 3:
	name = 'GRW70d5824'
	#this should be the file in AB magnitudes
	stdfile= '/Users/Briana/Documents/cure/virusp1/scripts/standard_stars/grw70d5824_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=119
else: 
	sys.exti("ID does not exist")

#--------------------#
# STANDARD VARIABLES #
#--------------------#

fileprefix = 'Spesvp'
DMprefix = 'pesvp'

num_fibs = 245 #246 but one not being found in trace for some reason
R0=float64(15) # radius in arcsec arround fiber0 to calculate barycenter
R1=float64(20) # radius in arcsec arround fiber0 to calculate barycenter
RF=8 # radius in arcsec to use for PSF modeling
aperture=5 # psf of aperture of fiber
halfap=float(aperture/2.0)
fib_radius = 2.15 #fiber radius is arcsec

#find the filename numbers from sci_im_lis.lis
sci_im_list = open('sci_im_lis.lis').read().splitlines()
filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

print 'FLUX CALIB FOR '+str(name)
print 'Files found: '+str(filelis)

#--------------------------------------------#
# FIND OBSERVATION PARAMETERS FROM DITH FILE #
#--------------------------------------------#

#name of the dither file (should be the same for all objects)
dith_file = 'dith.txt'

shift_x, shift_y, see_lis, airmass_lis = loadtxt(dith_file, usecols = (4,5,6,8), unpack=True)
seeing = average(see_lis) # seeing in Fiber Coordinate units, presumably arcseconds
airmass = average(airmass_lis)

print 'Seeing: '+str(seeing)+', Airmass: '+str(airmass)

#------------------------#
# GET STANDARD STAR DATA #
#------------------------#

print '[1] GETTING STANDARD STAR DATA'

#it is using the file that is in AB magnitudes
sl1, sm1 = loadtxt(stdfile, usecols = (0,1), unpack=True)

#mask out values greater than 7000A and less than 3400A for the interpolate routine in airmass correction
mask_wl = ma.masked_outside(sl1, 3400, 7000)
sl = mask_wl.compressed()
sm = ma.masked_array(sm1, mask=mask_wl.mask).compressed()

sfnu=pow(10.0,(-0.4*(sm+48.59))) # flux denisty in ergs/s/cm2/Hz
c=299792458.0e10 # c in A/s
sflam=(c/(pow(sl,2)))*sfnu # flux density in ergs/s/cm2/A

#----------------------------------------#
# EXTINCT STANDARD STAR FLUX FOR AIRMASS #
#----------------------------------------#

print '[2] EXTINCT STANDARD STAR FLUX FOR AIRMASS'

elam, ecoeff = loadtxt(op.join(curebin_path,'../scripts/extinc_McD.dat'), usecols = (0,1), unpack=True)
f = interpolate.interp1d(x=elam*10.0, y=ecoeff)
ecoeffint = f(sl)
sflam=sflam*(10.0**(-0.4*ecoeffint*airmass))

#------------------------------------#
# EXTRACT 1D SPECTRUM FOR EACH FIBER #
#------------------------------------#

#subtractsky -p S -d pesvp0073.dist -f pesvp0073.fmod -J -mask-5577  -w 200 pesvp0073.fits
#apimage -b M4861 -p Hbeta -d ../../redux/ap2015/20150417/morn_mastertrace.dist -i IFUcen_VP2_27m.txt -c dith.txt

if fiberextract == True:
	print '[3] PERFORMING FIBER EXTRACTION'
	for i in range(len(filelis)):
		bash_com = op.join(curebin_path,'fiberextract -d '+str(DMprefix)+str(filelis[i])+'.dist -f '+str(DMprefix)+str(filelis[i])+'.fmod -c -r 1,246 -l ['+str(wl_lower)+','+str(wl_higher)+'] '+str(fileprefix)+str(filelis[i])+'.fits')
		#$CUREBIN/fiberextract -d pesvp0068.dist -f pesvp0068.fmod -c -r 1,246 -l [3500,5800] Spesvp0068.fits 
		subprocess.call(bash_com, shell=True)

print 'Done with Fiber Extraction [3]'

#---------------------------------#
# READ FIBER EXTRACTED FITS FILES #
#---------------------------------#

print '[4] READING FIBER EXTRACTED FITS FILES'

im1 = pyfits.open('Fe'+str(fileprefix)+str(filelis[0])+'.fits')
d1 =  im1[0].data
h1 = im1[0].header

im2 = pyfits.open('Fe'+str(fileprefix)+str(filelis[1])+'.fits')
d2 =  im2[0].data
h2 = im2[0].header

im3 = pyfits.open('Fe'+str(fileprefix)+str(filelis[2])+'.fits')
d3 =  im3[0].data
h3 = im3[0].header

im4 = pyfits.open('Fe'+str(fileprefix)+str(filelis[3])+'.fits')
d4 =  im4[0].data
h4 = im4[0].header

im5 = pyfits.open('Fe'+str(fileprefix)+str(filelis[4])+'.fits')
d5 =  im5[0].data
h5 = im5[0].header

im6 = pyfits.open('Fe'+str(fileprefix)+str(filelis[5])+'.fits')
d6 =  im6[0].data
h6 = im6[0].header

#check if it is red and if it is cut of the last few pixels because they contain nan values
if wl_higher > 6000:
	print 'Trimming Red data to get rid of Nans and Infs'
	d1 = vstack((d1[0:244,0:900],zeros(900)))
	d2 = vstack((d2[0:244,0:900],zeros(900)))
	d3 = vstack((d3[0:244,0:900],zeros(900)))
	d4 = vstack((d4[0:244,0:900],zeros(900)))
	d5 = vstack((d5[0:244,0:900],zeros(900)))
	d6 = vstack((d6[0:244,0:900],zeros(900)))

di1 = count_nonzero(isnan(d1)) + count_nonzero(isinf(d1))
di2 = count_nonzero(isnan(d2)) + count_nonzero(isinf(d2))
di3 = count_nonzero(isnan(d3)) + count_nonzero(isinf(d3))
di4 = count_nonzero(isnan(d4)) + count_nonzero(isinf(d4))
di5 = count_nonzero(isnan(d5)) + count_nonzero(isinf(d5))
di6 = count_nonzero(isnan(d6)) + count_nonzero(isinf(d6))
print 'Number of Invalids Found in dithers: '+str(di1)+', '+str(di2)+', '+str(di3)+', '+str(di4)+', '+str(di5)+', '+str(di6)

exptime = int(h1['EXPTIME'])

#find length of one spectrum in a fiber
ln = len(d1[0])

#--------------------------------#
# READ COORDINATES OF EACH FIBER #
#--------------------------------#

print '[5] READING IFU CEN FILE'

#using the IFU mapping file to find the positons of the fibers for each dither
fn, fx1, fy1 = loadtxt(op.join(curebin_path,'../config/IFUcen_vp2_27m.txt'), usecols = (0,1,2), skiprows=3, unpack=True)
fx2=fx1+shift_x[1]
fy2=fy1+shift_y[1]   #Check that I changed the dither pattern correctly!
fx3=fx2+shift_x[2]
fy3=fy2+shift_y[2]
fx4=fx3+shift_x[3]
fy4=fy3+shift_y[3]
fx5=fx4+shift_x[4]
fy5=fy4+shift_y[4]
fx6=fx5+shift_x[5]
fy6=fy5+shift_y[5]

if create_standSpec:

	#-------------------------------------#
	# FIND THE AVERAGE CENTER OF THE STAR #
	#-------------------------------------#

	print '[6] FINDING THE AVERAGE CENTER OF THE STAR'

	#making an alias for the original data files because these will be manipulated 
	aux1=d1[0:num_fibs] 
	aux2=d2[0:num_fibs] 
	aux3=d3[0:num_fibs] 
	aux4=d4[0:num_fibs] 
	aux5=d5[0:num_fibs] 
	aux6=d6[0:num_fibs] 

	fx1=fx1[0:num_fibs]
	fy1=fy1[0:num_fibs]
	fx2=fx2[0:num_fibs]
	fy2=fy2[0:num_fibs]
	fx3=fx3[0:num_fibs]
	fy3=fy3[0:num_fibs]
	fx4=fx4[0:num_fibs]
	fy4=fy4[0:num_fibs]
	fx5=fx5[0:num_fibs]
	fy5=fy5[0:num_fibs]
	fx6=fx6[0:num_fibs]
	fy6=fy6[0:num_fibs]

	#integrating all light in each fiber, creates arrays containing total flux in each fiber for each dither
	totflux1 = sum(aux1, axis = 1)
	totflux2 = sum(aux2, axis = 1)
	totflux3 = sum(aux3, axis = 1)
	totflux4 = sum(aux4, axis = 1)
	totflux5 = sum(aux5, axis = 1)
	totflux6 = sum(aux6, axis = 1)

	#define the threshold for saying there is enough flux for this to be a star fiber
	thresh=5*median([totflux1,totflux2,totflux3,totflux4,totflux5,totflux6])

	print 'Threshold flux: '+str(thresh)

	#This is finding where the position of the fiber is inside of the defined radius R0
	#but The flux in those fibers also has to be greater than the defined threshold value above
	#each sel array is an array of whether each value (fiber) is true or false
	sel1=((sqrt(abs(pow(subtract(fx1,fx1[fib0]),2))+abs(pow(subtract(fy1,fy1[fib0]),2))) < R0 ) & (totflux1 > thresh))
	sel2=((sqrt(abs(pow(subtract(fx2,fx1[fib0]),2))+abs(pow(subtract(fy2,fy1[fib0]),2))) < R0 ) & (totflux2 > thresh))
	sel3=((sqrt(abs(pow(subtract(fx3,fx1[fib0]),2))+abs(pow(subtract(fy3,fy1[fib0]),2))) < R0 ) & (totflux3 > thresh))
	sel4=((sqrt(abs(pow(subtract(fx4,fx1[fib0]),2))+abs(pow(subtract(fy4,fy1[fib0]),2))) < R0 ) & (totflux4 > thresh))
	sel5=((sqrt(abs(pow(subtract(fx5,fx1[fib0]),2))+abs(pow(subtract(fy5,fy1[fib0]),2))) < R0 ) & (totflux5 > thresh))
	sel6=((sqrt(abs(pow(subtract(fx6,fx1[fib0]),2))+abs(pow(subtract(fy6,fy1[fib0]),2))) < R0 ) & (totflux6 > thresh))

	#center of mass x value
	#an array[sel#] will return the values at the fibers that met the threshold criteria 
	wx1=sum(multiply(fx1[sel1],totflux1[sel1]))
	wx2=sum(multiply(fx2[sel2],totflux2[sel2]))
	wx3=sum(multiply(fx3[sel3],totflux3[sel3]))
	wx4=sum(multiply(fx4[sel4],totflux4[sel4]))
	wx5=sum(multiply(fx5[sel5],totflux5[sel5]))
	wx6=sum(multiply(fx6[sel6],totflux6[sel6]))

	#center of mass x value
	wy1=sum(multiply(fy1[sel1],totflux1[sel1])) 
	wy2=sum(multiply(fy2[sel2],totflux2[sel2]))
	wy3=sum(multiply(fy3[sel3],totflux3[sel3]))
	wy4=sum(multiply(fy4[sel4],totflux4[sel4]))
	wy5=sum(multiply(fy5[sel5],totflux5[sel5]))
	wy6=sum(multiply(fy6[sel6],totflux6[sel6]))

	#flux weighted position, output is the positon of the centroid [sx,sy] of star
	sx=sum([wx1,wx2,wx3,wx4,wx5,wx6])/sum([sum(totflux1[sel1]),sum(totflux2[sel2]),sum(totflux3[sel3]),sum(totflux4[sel4]),sum(totflux5[sel5]),sum(totflux6[sel6])])
	sy=sum([wy1,wy2,wy3,wy4,wy5,wy6])/sum([sum(totflux1[sel1]),sum(totflux2[sel2]),sum(totflux3[sel3]),sum(totflux4[sel4]),sum(totflux5[sel5]),sum(totflux6[sel6])])

	print 'position of centroid of star: '
	print '	x: '+str(sx)
	print '	y: '+str(sy)

	#------------------------------------------------#
	# ATMOSPHERIC DIFFERENTIAL REFRACTION CORRECTION #
	#------------------------------------------------#

	print '[7] FINDING CENTOID OFFSETS FOR ADR CORRECTION'

	#this will find the ADR correction by finding calculating the center of the star at every wavelength to calculate the shift in the star due to ADR effects
	#to attempt to somwhat 'smooth' the solution it will caluclate the center with a few wavelenths (pixels) on either side of the one being used. 

	#this is the factor that you want to smooth it by. 
	#will include +/- sfac pixels when calculating the center for that pixel (wavelength)
	sfac = 10
	szero = zeros((len(aux1),sfac))

	#This is appending #sfac zeros to the begginning and end of each dither array so the loop can handle the edges 
	aux_s1 = hstack((szero,aux1,szero))
	aux_s2 = hstack((szero,aux2,szero))
	aux_s3 = hstack((szero,aux3,szero))
	aux_s4 = hstack((szero,aux4,szero))
	aux_s5 = hstack((szero,aux5,szero))
	aux_s6 = hstack((szero,aux6,szero))

	#these will hold all of the x and y offsets and the thresholds for each wavelength. 
	#the thresholds are just stored as a check if you want to look at statisics 
	ADR_x = zeros(ln)
	ADR_y = zeros(ln)
	ADR_thresh = zeros(ln)


	for i in range(ln):

		#This starts the iteration at the first wavelength
		j = i + sfac

		#This creates an array for each dither of just one wavelength (pixel) +/- #sfac values
		waveflux1 = sum(aux_s1[:,j-sfac:j+sfac+1],axis=1)
		waveflux2 = sum(aux_s2[:,j-sfac:j+sfac+1],axis=1)
		waveflux3 = sum(aux_s3[:,j-sfac:j+sfac+1],axis=1)
		waveflux4 = sum(aux_s4[:,j-sfac:j+sfac+1],axis=1)
		waveflux5 = sum(aux_s5[:,j-sfac:j+sfac+1],axis=1)
		waveflux6 = sum(aux_s6[:,j-sfac:j+sfac+1],axis=1)

		#define the threshold for saying there is enough flux for this to be a star fiber
		thresh=5*median([waveflux1,waveflux2,waveflux3,waveflux4,waveflux5,waveflux6])
		ADR_thresh[i]=thresh

		#This is finding where the position of the fiber is inside of the defined radius R0
		#but The flux in those fibers also has to be greater than the defined threshold value above
		#each sel array is an array of whether each value (fiber) is true or false
		sel1=((sqrt(abs(pow(subtract(fx1,fx1[fib0]),2))+abs(pow(subtract(fy1,fy1[fib0]),2))) < R1 ) & (waveflux1 > thresh))
		sel2=((sqrt(abs(pow(subtract(fx2,fx1[fib0]),2))+abs(pow(subtract(fy2,fy1[fib0]),2))) < R1 ) & (waveflux2 > thresh))
		sel3=((sqrt(abs(pow(subtract(fx3,fx1[fib0]),2))+abs(pow(subtract(fy3,fy1[fib0]),2))) < R1 ) & (waveflux3 > thresh))
		sel4=((sqrt(abs(pow(subtract(fx4,fx1[fib0]),2))+abs(pow(subtract(fy4,fy1[fib0]),2))) < R1 ) & (waveflux4 > thresh))
		sel5=((sqrt(abs(pow(subtract(fx5,fx1[fib0]),2))+abs(pow(subtract(fy5,fy1[fib0]),2))) < R1 ) & (waveflux5 > thresh))
		sel6=((sqrt(abs(pow(subtract(fx6,fx1[fib0]),2))+abs(pow(subtract(fy6,fy1[fib0]),2))) < R1 ) & (waveflux6 > thresh))

		#center of mass x value
		#an array[sel#] will return the values at the fibers that met the threshold criteria 
		wx1=sum(multiply(fx1[sel1],waveflux1[sel1]))
		wx2=sum(multiply(fx2[sel2],waveflux2[sel2]))
		wx3=sum(multiply(fx3[sel3],waveflux3[sel3]))
		wx4=sum(multiply(fx4[sel4],waveflux4[sel4]))
		wx5=sum(multiply(fx5[sel5],waveflux5[sel5]))
		wx6=sum(multiply(fx6[sel6],waveflux6[sel6]))

		#center of mass x value
		wy1=sum(multiply(fy1[sel1],waveflux1[sel1])) 
		wy2=sum(multiply(fy2[sel2],waveflux2[sel2]))
		wy3=sum(multiply(fy3[sel3],waveflux3[sel3]))
		wy4=sum(multiply(fy4[sel4],waveflux4[sel4]))
		wy5=sum(multiply(fy5[sel5],waveflux5[sel5]))
		wy6=sum(multiply(fy6[sel6],waveflux6[sel6]))

		#flux weighted position, output is the positon of the centroid [sx,sy] of star
		sx=sum([wx1,wx2,wx3,wx4,wx5,wx6])/sum([sum(waveflux1[sel1]),sum(waveflux2[sel2]),sum(waveflux3[sel3]),sum(waveflux4[sel4]),sum(waveflux5[sel5]),sum(waveflux6[sel6])])
		sy=sum([wy1,wy2,wy3,wy4,wy5,wy6])/sum([sum(waveflux1[sel1]),sum(waveflux2[sel2]),sum(waveflux3[sel3]),sum(waveflux4[sel4]),sum(waveflux5[sel5]),sum(waveflux6[sel6])])

		ADR_x[i] = sx
		ADR_y[i] = sy

	#ADR x and y give you how the centroid of the star shift with wavelength
	#normalize the centriod values by the first value to get a general posisitonal shift as a function of wavelength
	#ADR x and y give you how the centroid of the star shift with wavelength
	print 'Avg ADR X: '+str(average(ADR_x))+'  Stddev: '+str(std(ADR_x))
	print 'Avg ADR Y: '+str(average(ADR_y))+'  Stddev: '+str(std(ADR_y))

	#normalize the centriod values by the first value to get a general posisitonal shift as a function of wavelength
	ADR_x_norm = subtract(ADR_x, ADR_x[0])
	ADR_y_norm = subtract(ADR_y, ADR_y[0])

	#plot the shift in the centriod as a function of wavelength 
	fig = plt.figure()
	ax1 = fig.add_subplot(111)

	#plot the shift in the centriod as a function of wavelength 
	COLOR='blue'
	NPOINTS = len(ADR_x)
	MAP='rainbow'

	# cm = plt.get_cmap(MAP)
	# ax1.set_prop_cycle([cm(1.*i/(NPOINTS-1)) for i in range(NPOINTS-1)])
	# for i in range(NPOINTS-1):
	#     ax1.plot(ADR_x_norm[i:i+2],ADR_y_norm[i:i+2])

	# circle1=plt.Circle((ADR_x_norm[0],ADR_y_norm[0]),.01,color='black')
	# fig.gca().add_artist(circle1)

	# fig.savefig('ADR_StarCenShift_'+str(name)+'.png')
	# plt.show()

	#------------------------------------#
	# MAKE 2D GAUSSIAN MODEL OF THE STAR #
	#------------------------------------#

	print '[8] BUILDING MODELS OF STANDARD STAR'

	Npts    = 1e5 # should be enough points
	tot_fibs = num_fibs*len(filelis) #total number of fibers = number of fibers times number of dithers

	fiber_x = hstack((fx1,fx1,fx2,fx4,fx5,fx6)) # list of fiber x postions from all dithers
	fiber_y = hstack((fy1,fy1,fy2,fy4,fy5,fy6)) # list of fiber y postions from all dithers

	wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
	start_wave = float(h1['CRVAL1']) #this is the starting wavelength value

	wave = add(arange(wave_sol,(wave_sol*ln)+wave_sol,wave_sol),start_wave)  # wavelength solution for each of the fibers (should be same across the board)
	fiber_flux   = vstack((aux1,aux2,aux3,aux4,aux5,aux6)) # This is where all of the fiber fluxes are kept ()
	F            = zeros((len(wave),)) # least squares method

	for j in range(len(wave)):

		if j%10 == 0:
			print 'Model '+str(j)+' of '+str(len(wave))

		P = multivariate_normal([ADR_x[j],ADR_y[j]],seeing * np.eye(2),(Npts,))
		#scatter(P[:,0],P[:,1],s=5,alpha=0.02,edgecolor='none')

		w = zeros((tot_fibs,1)) # weight for each fiber

		for i in xrange(tot_fibs):
			#finding the difference in position from the fiber center to each of the points in the star gaussian model
			sel  = ((fiber_x[i] - P[:,0])**2 + (fiber_y[i] - P[:,1])**2) < fib_radius**2
			w[i] = divide(sum(sel), Npts) 

		total_weight = sum(w)
		total_flux   = divide(sum(fiber_flux,axis=0),total_weight) # Gary's method

		F[j] = lstsq(w,fiber_flux[:,j:j+1])[0]

	save('fluxCalibSpec_'+str(name), F)
	save('WavelengthSolu_'+str(name), wave)

print 'Done with building ADR corrected spectrum [6][7][8][9]'

if create_sensFunc: 

	#-------------------------------------#
	# GET OBSERVED STANDARD STAR SPECTRUM #
	#-------------------------------------#

	#load in standard spectrum as derived above 
	stand_spec = load('fluxCalibSpec_'+str(name)+'.npy') #DN/A

	#get the wavelength values for this spectrum 
	wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
	start_wave = float(h1['CRVAL1']) #this is the starting wavelength value

	#divide by exposure time to get per second
	stand_spec_perpix = divide(stand_spec,exptime) #DN/s/pixel
	#convert standard spectrum to be per anstrom so it is invariant
	stand_spec = divide(stand_spec_perpix,wave_sol) #DN/s/A

	wave = add(arange(wave_sol,(wave_sol*ln)+wave_sol,wave_sol),start_wave)

	#-------------------------------#
	# CALCULATE SENISIVITY FUNCTION #
	#-------------------------------#

	print '[9] FINDING SENSITIVITY FUNCTION'

	#interpolate the know spectrum and create the spectrum over the wavelength array for the standar star spectrum
	f = interpolate.interp1d(x=sl, y=sflam)
	sflamint = f(wave)

	#smooth over both the standard and known spectrum before divison using a simple median smoothing 
	stand_spec_sm = signal.medfilt(stand_spec, kernel_size=101)
	sflamint_sm = signal.medfilt(sflamint, kernel_size=101)

	#divide to get the sensitivity funciton 
	sens_func = divide(stand_spec, sflamint) #DN*cm^2 / erg
	sens_func_sm = divide(stand_spec_sm, sflamint_sm) #smoothed version 
	
	save('SensitivityFunc_'+str(name), sens_func_sm)
	save('WavelengthSolu_'+str(name), wave)

	#------------#
	# PLOT STUFF #
	#------------#

	plt.subplot(3,1,1)
	plt.plot(wave, stand_spec, color='blue', label='Standard')
	plt.plot(wave, stand_spec_sm, color = 'red', label='Stand Sm')
	plt.title('Observed Spectrum')
	plt.ylabel('DN/s/pixel')
	#plt.legend()

	plt.subplot(3,1,2)
	plt.plot(wave, sflamint, color = 'red', label='Known Spec')
	plt.plot(wave, sflamint_sm, color = 'purple', label='Known Sm')
	plt.plot(wave, divide(stand_spec_sm,sens_func_sm), color = 'green', label='recovery')
	plt.title('Oke Spectrum')
	plt.ylabel('Flux (ergs/s/cm^2/A)')
	#plt.legend()

	plt.subplot(3,1,3)
	plt.plot(wave, sens_func, color = 'purple', label='Known Sm')
	plt.plot(wave, sens_func_sm, color = 'orange', label='Stand Sm')
	plt.title('Sensitivity Function (Observed/Oke)')
	plt.ylabel('DN*cm^2 / erg')
	plt.xlabel('Wavelength (A)')
	#plt.legend()
	#plt.show()

	plt.savefig('Flux_Calib_'+str(name)+'.png')



