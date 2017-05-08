import numpy as np
import string
import matplotlib.pyplot as plt
import pyfits
import subprocess
from pylab import *
from scipy import interpolate

#run this script inside of the data folder 

#************************#
# User Defined Variables #
#************************#

#define ID to perform flux calibration on
field = 'M82_F2'
id_stand = 1
Trim_setting = 'apr15_red'

#prefix for your sciene files 
prefix = 'vp'

#store the flux calibrated array as a new fits file 
save_fits = False

#choose one fiber to plot in the end (123 is middle fiber)
plot_fib = 125

#----------------------------#
# User defined standard star #
#----------------------------#

#import sensitivity function
if id_stand == 1:
	name_stand = 'Feige34'
elif id_stand == 2: 
	name_stand = 'Feige67'
elif id_stand == 3:
	name_stand = 'GRW70d5824'
elif id_stand == 4:
	name_stand = 'BD_75D325'
else: 
	print "ID does not exist"

#-------------------------#
# User defined trim image #
#-------------------------#

if Trim_setting == 'apr16_03':
	x_srt = 0 	#default 0
	x_end = 990 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr16_04':
	x_srt = 0 	#default 0
	x_end = 920 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'apr16_05':
	x_srt = 0 	#default 0
	x_end = 985 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr16_06':
	x_srt = 0 	#default 0
	x_end = 910 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'apr16_07':
	x_srt = 0 	#default 0
	x_end = 910 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'mar17_30':
	x_srt = 12 	#default 0
	x_end = 912 #default 1024
	y_srt = 0   #default 0
	y_end = 243 #default 245
elif Trim_setting == 'mar17_31':
	x_srt = 34 	#default 0
	x_end = 988 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'feb17':
	x_srt = 35 	#default 0
	x_end = 988 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'mar16_red':
	x_srt = 20 	#default 0
	x_end = 920 #default 1024
	y_srt = 0   #default 0
	y_end = 243 #default 245
elif Trim_setting == 'mar16_blue':
	x_srt = 35 	#default 0
	x_end = 985 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr15_red':
	x_srt = 11 	#default 0
	x_end = 958 #default 1024
	y_srt = 0   #default 0
	y_end = 243 #default 245
elif Trim_setting == 'apr15_blue':
	x_srt = 34 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
else:
	print 'No trim setting applied'
	x_srt = 0 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245

#+++++++++++++++++++ end user defined variables +++++++++++++++++++++#

#----------#
# Get Data #
#----------#

#find the filename numbers from sci_im_lis.lis
sci_im_list = open(field+'/sci_im_lis.lis').read().splitlines()
filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

#--------------------------------------------#
# FIND OBSERVATION PARAMETERS FROM DITH FILE #
#--------------------------------------------#

see_lis, airmass_lis = np.loadtxt(field+'/dith.txt', usecols = (6,8), unpack=True)
seeing  = np.average(see_lis) # seeing in Fiber Coordinate units, presumably arcseconds
airmass = np.average(airmass_lis)

print 'Seeing: '+str(seeing)+', Airmass: '+str(airmass)

#--------------------------#
# LOAD SENSITIVTY FUNCTION #
#--------------------------#

sens_func_file = 'SensitivityFunc_'+name_stand+'.npy'
#load in the sensitivity function
sens_func = np.load(name_stand+'/'+sens_func_file)

#-----------------------------------#
# LOAD STANDARD WAVELENGTH SOLUTION #
#-----------------------------------#

wave_file = 'WavelengthSolu_'+name_stand+'.npy'
#load in the wavelength solution for standard star
wave_stand = np.load(name_stand+'/'+wave_file)

#---------------------------------#
# READ FIBER EXTRACTED FITS FILES #
#---------------------------------#

#import the fiber extracted files for the object 
im1 = pyfits.open(field+'/'+'FeSpes'+str(prefix)+str(filelis[0])+'.fits')
d1 =  im1[0].data
h1 = im1[0].header

im2 = pyfits.open(field+'/'+'FeSpes'+str(prefix)+str(filelis[1])+'.fits')
d2 =  im2[0].data
h2 = im2[0].header

im3 = pyfits.open(field+'/'+'FeSpes'+str(prefix)+str(filelis[2])+'.fits')
d3 =  im3[0].data
h3 = im3[0].header

im4 = pyfits.open(field+'/'+'FeSpes'+str(prefix)+str(filelis[3])+'.fits')
d4 =  im4[0].data
h4 = im4[0].header

im5 = pyfits.open(field+'/'+'FeSpes'+str(prefix)+str(filelis[4])+'.fits')
d5 =  im5[0].data
h5 = im5[0].header

im6 = pyfits.open(field+'/'+'FeSpes'+str(prefix)+str(filelis[5])+'.fits')
d6 =  im6[0].data
h6 = im6[0].header

#find length of one spectrum in a fiber before it is trimmed 
ln_orig = len(d1[0])

#--------------------------------------#
# BUILD WAVELENGTH SOLUTION FOR OBJECT #
#--------------------------------------#

wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
start_wave = float(h1['CRVAL1']) #this is the starting wavelength value
wave_orig = np.add(np.arange(wave_sol,(wave_sol*ln_orig)+(0.5*wave_sol),wave_sol),start_wave)  # wavelength solution for each of the fibers (should be same across the board)
wave = wave_orig[x_srt:x_end] 	#trim the wavelength solution to make the trimmed image

np.save(field+'/WavelengthSolu_'+str(field), wave)

print 'Wave Sol for Data:     '+str(np.amax(wave))+' ; '+str(np.amin(wave))
print 'Wave Sol for Standard: '+str(np.amax(wave_stand))+' ; '+str(np.amin(wave_stand))

#need to check if the wavelength solution of the data is outside the bounds of the sensitivity function
if np.amax(wave) > np.amax(wave_stand):
	print 'Data outside bounds of sensitivity function: Data Max WL: '+str(np.amax(wave))+' ; '+'Stand Max WL: '+str(np.amax(wave_stand))
	ang_diff = np.amax(wave) - np.amax(wave_stand) #This is the difference in the solutions in angstroms
	pix_diff = int(round(ang_diff / wave_sol)) #divide the angstroms to get pixel difference (round to make it an int)
	x_end = (x_end - (pix_diff+1)) #change the choosen x_end pixel to cut out data that can't be flux calibrated 
	wave = wave_orig[x_srt:x_end] 	#trim the wavelength solution to make the trimmed image
	print 'Trimming X by and extra '+str(pix_diff+1)+' pixels at the end'
	print 'New Wave Sol for Data:     '+str(np.amax(wave))+' ; '+str(np.amin(wave))

if np.amin(wave) < np.amin(wave_stand):
	print 'Data outside bounds of sensitivity function: Data Min WL: '+str(np.amin(wave))+' ; '+'Stand Min WL: '+str(np.amin(wave_stand))
	ang_diff = np.amin(wave_stand) - np.amin(wave) #This is the difference in the solutions in angstroms
	pix_diff = int(round(ang_diff / wave_sol)) #divide the angstroms to get pixel difference (round to make it an int)
	x_srt = (x_srt + (pix_diff+1)) #change the choosen x_end pixel to cut out data that can't be flux calibrated 
	wave = wave_orig[x_srt:x_end] 	#trim the wavelength solution to make the trimmed image
	print 'Trimming X by and extra '+str(pix_diff+1)+' pixels in the beginning'
	print 'New Wave Sol for Data:     '+str(np.amax(wave))+' ; '+str(np.amin(wave))

#interpolate over the sensitivity funciton based on its wavelength solution
f_stand = interpolate.interp1d(x=wave_stand, y=sens_func)
#build a sensitivity function based on the wavelength solution for the object spectra 
sens_func_interp = f_stand(wave)

#-----------#
# TRIM DATA #
#-----------#

extra_row = np.zeros(x_end - x_srt)

if y_end == 245:
	d1 = d1[y_srt:y_end,x_srt:x_end]
	d2 = d2[y_srt:y_end,x_srt:x_end]
	d3 = d3[y_srt:y_end,x_srt:x_end]
	d4 = d4[y_srt:y_end,x_srt:x_end]
	d5 = d5[y_srt:y_end,x_srt:x_end]
	d6 = d6[y_srt:y_end,x_srt:x_end]
elif y_end == 244:
	d1 = np.vstack((d1[y_srt:y_end,x_srt:x_end],extra_row))
	d2 = np.vstack((d2[y_srt:y_end,x_srt:x_end],extra_row))
	d3 = np.vstack((d3[y_srt:y_end,x_srt:x_end],extra_row))
	d4 = np.vstack((d4[y_srt:y_end,x_srt:x_end],extra_row))
	d5 = np.vstack((d5[y_srt:y_end,x_srt:x_end],extra_row))
	d6 = np.vstack((d6[y_srt:y_end,x_srt:x_end],extra_row))
elif y_end == 243:
	d1 = np.vstack((d1[y_srt:y_end,x_srt:x_end],np.vstack((extra_row,extra_row))))
	d2 = np.vstack((d2[y_srt:y_end,x_srt:x_end],np.vstack((extra_row,extra_row))))
	d3 = np.vstack((d3[y_srt:y_end,x_srt:x_end],np.vstack((extra_row,extra_row))))
	d4 = np.vstack((d4[y_srt:y_end,x_srt:x_end],np.vstack((extra_row,extra_row))))
	d5 = np.vstack((d5[y_srt:y_end,x_srt:x_end],np.vstack((extra_row,extra_row))))
	d6 = np.vstack((d6[y_srt:y_end,x_srt:x_end],np.vstack((extra_row,extra_row))))
else:
	sys.exit("Excluding too many fibers. This script not equipt to handel this")

di1 = np.count_nonzero(np.isnan(d1)) + np.count_nonzero(np.isinf(d1))
di2 = np.count_nonzero(np.isnan(d2)) + np.count_nonzero(np.isinf(d2))
di3 = np.count_nonzero(np.isnan(d3)) + np.count_nonzero(np.isinf(d3))
di4 = np.count_nonzero(np.isnan(d4)) + np.count_nonzero(np.isinf(d4))
di5 = np.count_nonzero(np.isnan(d5)) + np.count_nonzero(np.isinf(d5))
di6 = np.count_nonzero(np.isnan(d6)) + np.count_nonzero(np.isinf(d6))
print 'Number of Invalids Found in dithers: '+str(di1)+', '+str(di2)+', '+str(di3)+', '+str(di4)+', '+str(di5)+', '+str(di6)

tot_nonzeros = di1+di2+di3+di4+di5+di6
if tot_nonzeros > 0:
	sys.exit(str(tot_nonzeros)+" Nans or Infs found!")

dith_lis = [d1,d2,d3,d4,d5,d6]
head_lis = [h1,h2,h3,h4,h5,h6]
num_fibs = np.shape(d1)[0]
exptime = int(h1['EXPTIME']) #this is the exposure time of an image

#find length of one spectrum in a fiber after it is trimmed
ln = len(d1[0])


#-----------------------------------#
# GET EXTINCTION CURVE FOR MCDONALD #
#-----------------------------------#

#load in the extinction curve for McDonald and use airmass derived from the dither file 
elam, ecoeff = np.loadtxt('/Users/Briana/Documents/cure/virusp1/scripts/extinc_McD.dat', usecols = (0,1), unpack=True)
f = interpolate.interp1d(x=elam*10.0, y=ecoeff)
ecoeffint = f(wave)

#-------------------------------------------#
# PERFORM FLUX CALIBRATION ON EACH SPECTRUM #
#-------------------------------------------#

#np.divide the object spec by the sensitivity function to flux calibrate it 
#have to iterate through each spectra in each dither
store_orig = []
store_fcalib = []
for i in range(len(filelis)):
	dith = dith_lis[i]
	for f in range(num_fibs):

		spec = dith[f] #DN/pixel

		#extinct standard star for airmass sense the sensitivity function is not dependent on airmass
		spec = spec*(10.0**(-0.4*ecoeffint*airmass))

		#np.divide spectrum by exposure time to get per sec
		spec_perpix = np.divide(spec,exptime) #DN/s/pixel
		#np.divide spectrum by wavelength per pixel to get per A b/c this is invariant 
		spec_perA = np.divide(spec_perpix,wave_sol) #DN/s/A

		#np.divide by the sensitivity function to get the flux calibrated spectrum 
		fluxcal_spec = np.divide(spec_perA,sens_func_interp)

		#replace the spectrum in the array with the flux calibrated one
		dith[f] = fluxcal_spec

		#this just stores the flux calibrated and original spectrum of the chosen fiber for plotting
		if f == plot_fib:
			store_orig.append(spec_perA)
			store_fcalib.append(fluxcal_spec)

	if save_fits:
		#store the flux calibrated array as a new fits file 
		#hdu = pyfits.PrimaryHDU(dith)
		pyfits.writeto(field+'/'+'FcalFeSpes'+str(prefix)+str(filelis[i])+'.fits', dith, header = head_lis[i],clobber=True)
		print 'Created file: '+field+'/'+'FcalFeSpes'+str(prefix)+str(filelis[i])+'.fits'

		err_im  = pyfits.open(field+'/e.Spes'+str(prefix)+str(filelis[i])+'.fits')
		err_dat = err_im[0].data
		err_hdr = err_im[0].header
		pyfits.writeto(field+'/'+'e.FcalFeSpes'+str(prefix)+str(filelis[i])+'.fits', err_dat, header = err_hdr,clobber=True)
		print 'Copied file: '+field+'/'+'e.Spes'+str(prefix)+str(filelis[i])+'.fits to: '+field+'/'+'e.FcalFeSpes'+str(prefix)+str(filelis[i])+'.fits'

plt.subplot(2,1,1)
plt.plot(wave, store_orig[0], color = 'red', label='1')
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
plt.savefig(field+'/Flux_Calib_'+str(field)+'.png')



