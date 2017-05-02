import numpy as np
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

#This should be run inside of the standard star folder (data/stand_star)

#************************#
# User Defined Variables #
#************************#

#define ID to perform flux calibration on
id = 4
fiberextract = False
create_standSpec = True
create_sensFunc = True

#Set to None if you do not want to trim
Trim_setting = 'feb17'

#define wavelength range (in angstroms)
#wl_lower  = 4700
#wl_higher = 7200
wl_lower  = 3500
wl_higher = 5900

#set up path to your cure bin
#extinc_McD.dat needs to be saved in the scripts directory in your CURE folder
curebin_path  = '/Users/bindahl/Documents/cure/virusp1/bin'
vp_redux_path = '/Users/Briana/Documents/Grad_School/VIRUS_P/VP_reduction'

#Data name prefix
prefix = 'vp'

#--------------------------------------#
# User defined standard star variables #
#--------------------------------------#

#define all the files to use and correspoinding error files
if id == 1:
	name = 'Feige34'
	#this should be the file in AB magnitudes
	stdfile= 'feige_34_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=123
	#define all the files to use and correspoinding error files
elif id == 2:
	name = 'Feige67'
	#this should be the file in AB magnitudes
	stdfile= 'feige_67_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=123
elif id == 3:
	name = 'GRW70d5824'
	#this should be the file in AB magnitudes
	stdfile= 'grw70d5824_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=119 #for apr16 #150 for mar17
elif id == 4:
	name = 'BD_75d325'
	#this should be the file in AB magnitudes
	stdfile= 'bd75d325_oke.txt'
	#y value of pixel you think is the center of the star in central dither (4th dither?)
	fib0=137 #for feb17 #151 
else: 
	sys.exti("ID does not exist")

#-------------------------#
# User defined trim image #
#-------------------------#

if Trim_setting == 'mar17_31':
	x_srt = 0 	#default 0
	x_end = 980 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'mar17_30':
	x_srt = 15 	#default 0
	x_end = 900 #default 1024
	y_srt = 0   #default 0
	y_end = 243 #default 245
elif Trim_setting == 'feb17':
	x_srt = 0 	#default 0
	x_end = 985 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr16_03':
	x_srt = 0 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr16_04':
	x_srt = 0 	#default 0
	x_end = 920 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'apr16_05':
	x_srt = 0 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'apr16_06':
	x_srt = 0 	#default 0
	x_end = 913 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'apr16_07':
	x_srt = 0 	#default 0
	x_end = 912 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
elif Trim_setting == 'mar16_red':
	x_srt = 25 	#default 0
	x_end = 915 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'mar16_blue':
	x_srt = 35 	#default 0
	x_end = 980 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr15_red':
	x_srt = 14 	#default 0
	x_end = 850 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
elif Trim_setting == 'apr15_blue':
	x_srt = 35 	#default 0
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

#--------------------#
# STANDARD VARIABLES #
#--------------------#

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

shift_x, shift_y, see_lis, airmass_lis = np.loadtxt(dith_file, usecols = (4,5,6,8), unpack=True)
seeing = np.average(see_lis) # seeing in Fiber Coordinate units, presumably arcseconds
airmass = np.average(airmass_lis)

print 'Seeing: '+str(seeing)+', Airmass: '+str(airmass)

#------------------------#
# GET STANDARD STAR DATA #
#------------------------#

print '[1] GETTING STANDARD STAR DATA'

#it is using the file that is in AB magnitudes
sl1, sm1 = np.loadtxt(op.join(vp_redux_path,'standard_stars',stdfile), usecols = (0,1), unpack=True)

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

elam, ecoeff = np.loadtxt(op.join(vp_redux_path,'extinc_McD.dat'), usecols = (0,1), unpack=True)
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
		bash_com = op.join(curebin_path,'fiberextract -d pes'+str(prefix)+str(filelis[i])+'.dist -f pes'+str(prefix)+str(filelis[i])+'.fmod -c -r 1,246 -l ['+str(wl_lower)+','+str(wl_higher)+'] Spes'+str(prefix)+str(filelis[i])+'.fits')
		#$CUREBIN/fiberextract -d pesvp0068.dist -f pesvp0068.fmod -c -r 1,246 -l [3500,5800] Spesvp0068.fits 
		subprocess.call(bash_com, shell=True)

print 'Done with Fiber Extraction [3]'

#---------------------------------#
# READ FIBER EXTRACTED FITS FILES #
#---------------------------------#

print '[4] READING FIBER EXTRACTED FITS FILES'

im1 = pyfits.open('FeSpes'+str(prefix)+str(filelis[0])+'.fits')
d1 =  im1[0].data
h1 = im1[0].header

im2 = pyfits.open('FeSpes'+str(prefix)+str(filelis[1])+'.fits')
d2 =  im2[0].data
h2 = im2[0].header

im3 = pyfits.open('FeSpes'+str(prefix)+str(filelis[2])+'.fits')
d3 =  im3[0].data
h3 = im3[0].header

im4 = pyfits.open('FeSpes'+str(prefix)+str(filelis[3])+'.fits')
d4 =  im4[0].data
h4 = im4[0].header

im5 = pyfits.open('FeSpes'+str(prefix)+str(filelis[4])+'.fits')
d5 =  im5[0].data
h5 = im5[0].header

im6 = pyfits.open('FeSpes'+str(prefix)+str(filelis[5])+'.fits')
d6 =  im6[0].data
h6 = im6[0].header

#check if it is red and if it is cut of the last few pixels because they contain nan values
# if wl_higher > 6000:
# 	print 'Trimming Red data to get rid of Nans and Infs'
# 	d1 = np.vstack((d1[0:244,0:900],np.zeros(900)))
# 	d2 = np.vstack((d2[0:244,0:900],np.zeros(900)))
# 	d3 = np.vstack((d3[0:244,0:900],np.zeros(900)))
# 	d4 = np.vstack((d4[0:244,0:900],np.zeros(900)))
# 	d5 = np.vstack((d5[0:244,0:900],np.zeros(900)))
# 	d6 = np.vstack((d6[0:244,0:900],np.zeros(900)))

#find length of one spectrum in a fiber before it is trimmed 
ln_orig = len(d1[0])

#trim data based on user settings to get ride of NaNs and Infs 
d1 = d1[y_srt:y_end,x_srt:x_end]
d2 = d2[y_srt:y_end,x_srt:x_end]
d3 = d3[y_srt:y_end,x_srt:x_end]
d4 = d4[y_srt:y_end,x_srt:x_end]
d5 = d5[y_srt:y_end,x_srt:x_end]
d6 = d6[y_srt:y_end,x_srt:x_end]

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

exptime = int(h1['EXPTIME'])

#find length of one spectrum in a fiber after it is trimmed
ln = len(d1[0])

print "LEN: "+str(ln)

#--------------------------------#
# READ COORDINATES OF EACH FIBER #
#--------------------------------#

print '[5] READING IFU CEN FILE'

#using the IFU mapping file to find the positons of the fibers for each dither
fn, fx1, fy1 = np.loadtxt(op.join(vp_redux_path,'IFUcen_vp2_27m.txt'), usecols = (0,1,2), skiprows=3, unpack=True)
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
	aux1=d1[0:y_end] 
	aux2=d2[0:y_end] 
	aux3=d3[0:y_end] 
	aux4=d4[0:y_end] 
	aux5=d5[0:y_end] 
	aux6=d6[0:y_end] 

	fx1=fx1[0:y_end]
	fy1=fy1[0:y_end]
	fx2=fx2[0:y_end]
	fy2=fy2[0:y_end]
	fx3=fx3[0:y_end]
	fy3=fy3[0:y_end]
	fx4=fx4[0:y_end]
	fy4=fy4[0:y_end]
	fx5=fx5[0:y_end]
	fy5=fy5[0:y_end]
	fx6=fx6[0:y_end]
	fy6=fy6[0:y_end]

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
	sel1=((np.sqrt(abs(pow(np.subtract(fx1,fx1[fib0]),2))+abs(pow(np.subtract(fy1,fy1[fib0]),2))) < R0 ) & (totflux1 > thresh))
	sel2=((np.sqrt(abs(pow(np.subtract(fx2,fx1[fib0]),2))+abs(pow(np.subtract(fy2,fy1[fib0]),2))) < R0 ) & (totflux2 > thresh))
	sel3=((np.sqrt(abs(pow(np.subtract(fx3,fx1[fib0]),2))+abs(pow(np.subtract(fy3,fy1[fib0]),2))) < R0 ) & (totflux3 > thresh))
	sel4=((np.sqrt(abs(pow(np.subtract(fx4,fx1[fib0]),2))+abs(pow(np.subtract(fy4,fy1[fib0]),2))) < R0 ) & (totflux4 > thresh))
	sel5=((np.sqrt(abs(pow(np.subtract(fx5,fx1[fib0]),2))+abs(pow(np.subtract(fy5,fy1[fib0]),2))) < R0 ) & (totflux5 > thresh))
	sel6=((np.sqrt(abs(pow(np.subtract(fx6,fx1[fib0]),2))+abs(pow(np.subtract(fy6,fy1[fib0]),2))) < R0 ) & (totflux6 > thresh))

	#center of mass x value
	#an array[sel#] will return the values at the fibers that met the threshold criteria 
	wx1=sum(np.multiply(fx1[sel1],totflux1[sel1]))
	wx2=sum(np.multiply(fx2[sel2],totflux2[sel2]))
	wx3=sum(np.multiply(fx3[sel3],totflux3[sel3]))
	wx4=sum(np.multiply(fx4[sel4],totflux4[sel4]))
	wx5=sum(np.multiply(fx5[sel5],totflux5[sel5]))
	wx6=sum(np.multiply(fx6[sel6],totflux6[sel6]))

	#center of mass x value
	wy1=sum(np.multiply(fy1[sel1],totflux1[sel1])) 
	wy2=sum(np.multiply(fy2[sel2],totflux2[sel2]))
	wy3=sum(np.multiply(fy3[sel3],totflux3[sel3]))
	wy4=sum(np.multiply(fy4[sel4],totflux4[sel4]))
	wy5=sum(np.multiply(fy5[sel5],totflux5[sel5]))
	wy6=sum(np.multiply(fy6[sel6],totflux6[sel6]))

	#flux weighted position, output is the positon of the centroid [sx,sy] of star
	sx=sum([wx1,wx2,wx3,wx4,wx5,wx6])/sum([sum(totflux1[sel1]),sum(totflux2[sel2]),sum(totflux3[sel3]),sum(totflux4[sel4]),sum(totflux5[sel5]),sum(totflux6[sel6])])
	sy=sum([wy1,wy2,wy3,wy4,wy5,wy6])/sum([sum(totflux1[sel1]),sum(totflux2[sel2]),sum(totflux3[sel3]),sum(totflux4[sel4]),sum(totflux5[sel5]),sum(totflux6[sel6])])

	if np.isnan(sx) or np.isnan(sy):
		sys.exit('CAN\'T FIND THE STAR: please check the central fiber defined for your standard star (fib0)')

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
	szero = np.zeros((len(aux1),sfac))

	#This is appending #sfac zeros to the begginning and end of each dither array so the loop can handle the edges 
	aux_s1 = np.hstack((szero,aux1,szero))
	aux_s2 = np.hstack((szero,aux2,szero))
	aux_s3 = np.hstack((szero,aux3,szero))
	aux_s4 = np.hstack((szero,aux4,szero))
	aux_s5 = np.hstack((szero,aux5,szero))
	aux_s6 = np.hstack((szero,aux6,szero))

	#these will hold all of the x and y offsets and the thresholds for each wavelength. 
	#the thresholds are just stored as a check if you want to look at statisics 
	ADR_x = np.zeros(ln)
	ADR_y = np.zeros(ln)
	ADR_thresh = np.zeros(ln)


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
		thresh=5*np.median([waveflux1,waveflux2,waveflux3,waveflux4,waveflux5,waveflux6])
		ADR_thresh[i]=thresh

		#This is finding where the position of the fiber is inside of the defined radius R0
		#but The flux in those fibers also has to be greater than the defined threshold value above
		#each sel array is an array of whether each value (fiber) is true or false
		sel1=((np.sqrt(abs(pow(np.subtract(fx1,fx1[fib0]),2))+abs(pow(np.subtract(fy1,fy1[fib0]),2))) < R1 ) & (waveflux1 > thresh))
		sel2=((np.sqrt(abs(pow(np.subtract(fx2,fx1[fib0]),2))+abs(pow(np.subtract(fy2,fy1[fib0]),2))) < R1 ) & (waveflux2 > thresh))
		sel3=((np.sqrt(abs(pow(np.subtract(fx3,fx1[fib0]),2))+abs(pow(np.subtract(fy3,fy1[fib0]),2))) < R1 ) & (waveflux3 > thresh))
		sel4=((np.sqrt(abs(pow(np.subtract(fx4,fx1[fib0]),2))+abs(pow(np.subtract(fy4,fy1[fib0]),2))) < R1 ) & (waveflux4 > thresh))
		sel5=((np.sqrt(abs(pow(np.subtract(fx5,fx1[fib0]),2))+abs(pow(np.subtract(fy5,fy1[fib0]),2))) < R1 ) & (waveflux5 > thresh))
		sel6=((np.sqrt(abs(pow(np.subtract(fx6,fx1[fib0]),2))+abs(pow(np.subtract(fy6,fy1[fib0]),2))) < R1 ) & (waveflux6 > thresh))

		#center of mass x value
		#an array[sel#] will return the values at the fibers that met the threshold criteria 
		wx1=sum(np.multiply(fx1[sel1],waveflux1[sel1]))
		wx2=sum(np.multiply(fx2[sel2],waveflux2[sel2]))
		wx3=sum(np.multiply(fx3[sel3],waveflux3[sel3]))
		wx4=sum(np.multiply(fx4[sel4],waveflux4[sel4]))
		wx5=sum(np.multiply(fx5[sel5],waveflux5[sel5]))
		wx6=sum(np.multiply(fx6[sel6],waveflux6[sel6]))

		#center of mass x value
		wy1=sum(np.multiply(fy1[sel1],waveflux1[sel1])) 
		wy2=sum(np.multiply(fy2[sel2],waveflux2[sel2]))
		wy3=sum(np.multiply(fy3[sel3],waveflux3[sel3]))
		wy4=sum(np.multiply(fy4[sel4],waveflux4[sel4]))
		wy5=sum(np.multiply(fy5[sel5],waveflux5[sel5]))
		wy6=sum(np.multiply(fy6[sel6],waveflux6[sel6]))

		#flux weighted position, output is the positon of the centroid [sx,sy] of star
		sx=sum([wx1,wx2,wx3,wx4,wx5,wx6])/sum([sum(waveflux1[sel1]),sum(waveflux2[sel2]),sum(waveflux3[sel3]),sum(waveflux4[sel4]),sum(waveflux5[sel5]),sum(waveflux6[sel6])])
		sy=sum([wy1,wy2,wy3,wy4,wy5,wy6])/sum([sum(waveflux1[sel1]),sum(waveflux2[sel2]),sum(waveflux3[sel3]),sum(waveflux4[sel4]),sum(waveflux5[sel5]),sum(waveflux6[sel6])])

		ADR_x[i] = sx
		ADR_y[i] = sy

	#ADR x and y give you how the centroid of the star shift with wavelength
	#normalize the centriod values by the first value to get a general posisitonal shift as a function of wavelength
	#ADR x and y give you how the centroid of the star shift with wavelength
	print 'Avg ADR X: '+str(np.average(ADR_x))+'  Stddev: '+str(np.std(ADR_x))
	print 'Avg ADR Y: '+str(np.average(ADR_y))+'  Stddev: '+str(np.std(ADR_y))

	#normalize the centriod values by the first value to get a general posisitonal shift as a function of wavelength
	ADR_x_norm = np.subtract(ADR_x, ADR_x[0])
	ADR_y_norm = np.subtract(ADR_y, ADR_y[0])

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
	tot_fibs = y_end*len(filelis) #total number of fibers = number of fibers times number of dithers

	fiber_x = np.hstack((fx1,fx1,fx2,fx4,fx5,fx6)) # list of fiber x postions from all dithers
	fiber_y = np.hstack((fy1,fy1,fy2,fy4,fy5,fy6)) # list of fiber y postions from all dithers

	wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
	start_wave = float(h1['CRVAL1']) #this is the starting wavelength value

	wave_orig = np.add(np.arange(wave_sol,(wave_sol*ln_orig)+(0.5*wave_sol),wave_sol),start_wave)  # wavelength solution for each of the fibers (should be same across the board)
	wave = wave_orig[x_srt:x_end]	#trim the wavelength solution to make the trimmed image

	fiber_flux   = np.vstack((aux1,aux2,aux3,aux4,aux5,aux6)) # This is where all of the fiber fluxes are kept ()
	F            = np.zeros((len(wave),)) # least squares method

	for j in range(len(wave)):

		if j%10 == 0:
			print 'Model '+str(j)+' of '+str(len(wave))

		P = multivariate_normal([ADR_x[j],ADR_y[j]],seeing * np.eye(2),(Npts,))
		#scatter(P[:,0],P[:,1],s=5,alpha=0.02,edgecolor='none')

		w = np.zeros((tot_fibs,1)) # weight for each fiber

		for i in xrange(tot_fibs):
			#finding the difference in position from the fiber center to each of the points in the star gaussian model
			sel  = ((fiber_x[i] - P[:,0])**2 + (fiber_y[i] - P[:,1])**2) < fib_radius**2
			w[i] = np.divide(sum(sel), Npts) 

		total_weight = sum(w)
		total_flux   = np.divide(sum(fiber_flux,axis=0),total_weight) # Gary's method

		F[j] = lstsq(w,fiber_flux[:,j:j+1])[0] #more robust method that is used 

	np.save('fluxCalibSpec_'+str(name), F)
	np.save('WavelengthSolu_'+str(name), wave)

print 'Done with building ADR corrected spectrum [6][7][8][9]'

if create_sensFunc: 

	#-------------------------------------#
	# GET OBSERVED STANDARD STAR SPECTRUM #
	#-------------------------------------#

	#load in standard spectrum as derived above 
	stand_spec = np.load('fluxCalibSpec_'+str(name)+'.npy') #DN/A

	#get the wavelength values for this spectrum 
	wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
	start_wave = float(h1['CRVAL1']) #this is the starting wavelength value

	#np.divide by exposure time to get per second
	stand_spec_perpix = np.divide(stand_spec,exptime) #DN/s/pixel
	#convert standard spectrum to be per anstrom so it is invariant
	stand_spec = np.divide(stand_spec_perpix,wave_sol) #DN/s/A

	wave_orig = np.add(np.arange(wave_sol,(wave_sol*ln_orig)+(0.5*wave_sol),wave_sol),start_wave)
	wave = wave_orig[x_srt:x_end] 	#trim the wavelength solution to make the trimmed image

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

	#np.divide to get the sensitivity funciton 
	sens_func = np.divide(stand_spec, sflamint) #DN*cm^2 / erg
	sens_func_sm = np.divide(stand_spec_sm, sflamint_sm) #smoothed version 
	
	np.save('SensitivityFunc_'+str(name), sens_func_sm)
	np.save('WavelengthSolu_'+str(name), wave)

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
	plt.plot(wave, np.divide(stand_spec_sm,sens_func_sm), color = 'green', label='recovery')
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



