import numpy as np
import string
import matplotlib.pyplot as plt
from scipy import interpolate
import pyfits
import subprocess
import numpy.ma as ma
from pylab import *
import os.path as op
import sys

#************************#
# User Defined Variables #
#************************#
data_file 	= '/Users/Briana/Documents/Grad_School/VIRUS_P/data_cubes/CuFcalFeSpes_F1_B_vp0054.fits'
error_file  = '/Users/Briana/Documents/Grad_School/VIRUS_P/data_cubes/e.CuFcalFeSpes_F1_B_vp0054.fits'
wl_sol_file = '/Users/Briana/Documents/Grad_School/VIRUS_P/wave_sols/WavelengthSolu_M82_F1_B.npy'

wl 	= 4861 #in Angstroms 
z	= 0.000677

wl_range = 30 #in pixels away from peak
peak_vals = [10]

#+++++++++++++++++++ end user defined variables +++++++++++++++++++++#

#***********#
# Load Data #
#***********#

#load the data file 
im = pyfits.open(data_file)
dat = im[0].data
hdr = im[0].header

#load the error data file 
im_e = pyfits.open(error_file)
dat_e = im_e[0].data
hdr_e = im_e[0].header

#load the wavelength solution
wave_sol = np.load(wl_sol_file) #AA

#Find if the data file loaded is a data cube or a 2D spectrum by measuring the dimensions 
num_dims = len(dat.shape)
if num_dims == 3:
	is_cube = True
	print 'Found a data cube!'
else:
	is_cube = False
	print 'Found a 2D spectrum'

#**************#
# Find Spectra #
#**************#

if is_cube:
	#starter spectrum
	spec1 = dat[:,20,20]

	#make the values reasonable numbers
	spec_avg = np.median(spec1) #find median value of spectrum
	avg_mag  = int(math.log10(spec_avg)) #find the order of magnitude of the median value
	spec1 	 = np.divide(spec1,pow(10,avg_mag)) #divide the spectrum by 10^(avg_mag)
	print 'Divided Spectrum by 1e'+str(avg_mag)

	wl_idx = (np.abs(wave_sol-wl)).argmin() #finds index of value closest to peak wavelength value 
	spec1_trim    = spec1[wl_idx-wl_range, wl_idx+wl_range]
	wave_sol_trim = wave_sol[wl_idx-wl_range, wl_idx+wl_range]

	plt.subplot(1,1,1)
	plt.plot(wave_sol_trim, spec1_trim, color='blue', label='Standard')
	#plt.plot(wave, stand_spec_sm, color = 'red', label='Stand Sm')
	plt.title('Spectrum around '+str(wl)+' (A)')
	plt.ylabel('Flux (ergs/s/cm^2/A) * 1e'+str(avg_mag))
	plt.xlabel('Wavelength (A)')
	plt.show()

else:
	print 'you have a 2D spec'

	


