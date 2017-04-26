from numpy import *
import string
import matplotlib.pyplot as plt
import pyfits
import subprocess
from pylab import *
from scipy import interpolate

#run this script inside of the data folder 

#define ID to perform flux calibration on
name = 'M82_F8'
color = 'red'
id_stand = 3

save_fits = True

#choose one fiber to plot in the end (123 is middle fiber)
plot_fib = 125

#----------#
# Get Data #
#----------#

#import sensitivity function
if id_stand == 1:
	name_stand = 'Feige34'
	sens_func_file = 'SensitivityFunc_Feige34.npy'
elif id_stand == 2: 
	name_stand = 'Feige67'
	sens_func_file = 'SensitivityFunc_Feige67.npy'
elif id_stand == 3:
	name_stand = 'GRW70d5824'
	sens_func_file = 'SensitivityFunc_GRW70d5824.npy'

else: 
	print "ID does not exist"

#find the filename numbers from sci_im_lis.lis
sci_im_list = open(name+'/sci_im_lis.lis').read().splitlines()
filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

fileprefix = 'Spesvp'
DMprefix = 'pesvp'

#--------------------------------------------#
# FIND OBSERVATION PARAMETERS FROM DITH FILE #
#--------------------------------------------#

see_lis, airmass_lis = loadtxt(name+'/dith.txt', usecols = (6,8), unpack=True)
seeing = average(see_lis) # seeing in Fiber Coordinate units, presumably arcseconds
airmass = average(airmass_lis)

print 'Seeing: '+str(seeing)+', Airmass: '+str(airmass)

#--------------------------#
# LOAD SENSITIVTY FUNCTION #
#--------------------------#

#load in the sensitivity function
sens_func = load(name_stand+'/'+sens_func_file)

#---------------------------------#
# READ FIBER EXTRACTED FITS FILES #
#---------------------------------#

#import the fiber extracted files for the object 
im1 = pyfits.open(name+'/'+'Fe'+str(fileprefix)+str(filelis[0])+'.fits')
d1 =  im1[0].data
h1 = im1[0].header

im2 = pyfits.open(name+'/'+'Fe'+str(fileprefix)+str(filelis[1])+'.fits')
d2 =  im2[0].data
h2 = im2[0].header

im3 = pyfits.open(name+'/'+'Fe'+str(fileprefix)+str(filelis[2])+'.fits')
d3 =  im3[0].data
h3 = im3[0].header

im4 = pyfits.open(name+'/'+'Fe'+str(fileprefix)+str(filelis[3])+'.fits')
d4 =  im4[0].data
h4 = im4[0].header

im5 = pyfits.open(name+'/'+'Fe'+str(fileprefix)+str(filelis[4])+'.fits')
d5 =  im5[0].data
h5 = im5[0].header

im6 = pyfits.open(name+'/'+'Fe'+str(fileprefix)+str(filelis[5])+'.fits')
d6 =  im6[0].data
h6 = im6[0].header

#check if it is red and if it is cut of the last few pixels because they contain nan values
if color == 'red':
	print 'Trimming Red data to get rid of Nans and Infs'
	d1 = vstack((d1[0:244,0:900],zeros(900)))
	d2 = vstack((d2[0:244,0:900],zeros(900)))
	d3 = vstack((d3[0:244,0:900],zeros(900)))
	d4 = vstack((d4[0:244,0:900],zeros(900)))
	d5 = vstack((d5[0:244,0:900],zeros(900)))
	d6 = vstack((d6[0:244,0:900],zeros(900)))

dith_lis = [d1,d2,d3,d4,d5,d6]
head_lis = [h1,h2,h3,h4,h5,h6]
num_fibs = shape(d1)[0]
exptime = int(h1['EXPTIME']) #this is the exposure time of an image

#find length of one spectrum in a fiber
ln = len(d1[0])

wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
start_wave = float(h1['CRVAL1']) #this is the starting wavelength value
wave = add(arange(wave_sol,(wave_sol*ln)+wave_sol,wave_sol),start_wave)  # wavelength solution for each of the fibers (should be same across the board)

#-----------------------------------#
# GET EXTINCTION CURVE FOR MCDONALD #
#-----------------------------------#

#load in the extinction curve for McDonald and use airmass derived from the dither file 
elam, ecoeff = loadtxt('/Users/Briana/Documents/cure/virusp1/scripts/extinc_McD.dat', usecols = (0,1), unpack=True)
f = interpolate.interp1d(x=elam*10.0, y=ecoeff)
ecoeffint = f(wave)

#-------------------------------------------#
# PERFORM FLUX CALIBRATION ON EACH SPECTRUM #
#-------------------------------------------#

#divide the object spec by the sensitivity function to flux calibrate it 
#have to iterate through each spectra in each dither
store_orig = []
store_fcalib = []
for i in range(len(filelis)):
	dith = dith_lis[i]
	for f in range(num_fibs):

		spec = dith[f] #DN/pixel

		#extinct standard star for airmass sense the sensitivity function is not dependent on airmass
		spec = spec*(10.0**(-0.4*ecoeffint*airmass))

		#divide spectrum by exposure time to get per sec
		spec_perpix = divide(spec,exptime) #DN/s/pixel
		#divide spectrum by wavelength per pixel to get per A b/c this is invariant 
		spec_perA = divide(spec_perpix,wave_sol) #DN/s/A

		#divide by the sensitivity function to get the flux calibrated spectrum 
		fluxcal_spec = divide(spec_perA,sens_func)

		#replace the spectrum in the array with the flux calibrated one
		dith[f] = fluxcal_spec

		#this just stores the flux calibrated and original spectrum of the chosen fiber for plotting
		if f == plot_fib:
			store_orig.append(spec_perA)
			store_fcalib.append(fluxcal_spec)

	if save_fits:
		#store the flux calibrated array as a new fits file 
		#hdu = pyfits.PrimaryHDU(dith)
		pyfits.writeto(name+'/'+'FcalFe'+str(fileprefix)+str(filelis[i])+'.fits', dith, header = head_lis[i],clobber=True)
		print 'Created file: '+name+'/'+'FcalFe'+str(fileprefix)+str(filelis[i])+'.fits'

plt.subplot(2,1,1)
plt.plot(wave, store_orig[0], color='red', label='1')
plt.plot(wave, store_orig[1], color = 'orange', label='2')
plt.plot(wave, store_orig[2], color = 'green', label='2')
plt.plot(wave, store_orig[3], color = 'blue', label='2')
plt.plot(wave, store_orig[4], color = 'purple', label='2')
plt.plot(wave, store_orig[5], color = 'black', label='2')
plt.title('Original Spectrum: Fiber '+str(plot_fib))
plt.ylabel('DN/s/A')

plt.subplot(2,1,2)
plt.plot(wave, store_fcalib[0], color='red', label='1')
plt.plot(wave, store_fcalib[1], color = 'orange', label='2')
plt.plot(wave, store_fcalib[2], color = 'green', label='2')
plt.plot(wave, store_fcalib[3], color = 'blue', label='2')
plt.plot(wave, store_fcalib[4], color = 'purple', label='2')
plt.plot(wave, store_fcalib[5], color = 'black', label='2')
plt.title('Flux Calibrated Spectrum: Fiber '+str(plot_fib))
plt.ylabel('ergs/s/cm^2/A')

#plt.show()
plt.savefig(name+'/Flux_Calib_'+str(name)+'.png')



