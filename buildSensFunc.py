import numpy as np
import numpy.ma as ma
import string
import matplotlib.pyplot as plt
from pylab import *
from scipy import interpolate
from scipy import signal 
from scipy import stats as st
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import pyfits
import subprocess
import os.path as op
import sys
import glob
from datetime import datetime
from astropy.coordinates import Angle

#This should be run inside of the standard star folder (data/stand_star)

#************************#
# User Defined Variables #
#************************#

#define ID to perform flux calibration on
fiberextract = False
create_standSpec = True
create_sensFunc = True

#Set to None if you do not want to trim
Trim_setting = 'apr16_06'

#set up path to your cure bin
#extinc_McD.dat needs to be saved in the scripts directory in your CURE folder
curebin_path  = '/Users/bindahl/Documents/cure/virusp1/bin'
vp_redux_path = '/Users/Briana/Documents/Grad_School/VIRUS_P/VP_reduction'

#path to wrk file with object coordiantes 
coord_file = '/Volumes/Briana_mac3/M82_data/VP_M82_data/object_lis.txt'

#Data name prefix
prefix = 'vp'

#name of the dither file (should be the same for all objects)
dith_file = 'dith.txt'

#find the filename numbers from sci_im_lis.lis
sci_im_list = glob.glob('FeSp*')

#latitude of McDonald Observatory
#McD_lat = 30.6798
McD_lat = 33.1
#-------------------------#
# User defined trim image #
#-------------------------#

if Trim_setting == 'mar17_31':
	id = 4 #define standard star id to use (see defined below)
	wl_lower  = 3500 #define wavelength range (in angstroms)
	wl_higher = 5900 
	x_srt = 0 	#default 0
	x_end = 985 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
	weath = 'mar2017' #define the McD weather file to use 
elif Trim_setting == 'mar17_30':
	id = 4
	wl_lower  = 4700 
	wl_higher = 7200
	x_srt = 12 	#default 0
	x_end = 912 #default 1024
	y_srt = 0   #default 0
	y_end = 243 #default 245
	weath = 'mar2017'
elif Trim_setting == 'feb17':
	id = 4
	wl_lower  = 3500 
	wl_higher = 5900 
	x_srt = 0 	#default 0
	x_end = 985 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
	weath = 'feb2017'
elif Trim_setting == 'apr16_03':
	id = 1
	wl_lower  = 3500 
	wl_higher = 5900 
	x_srt = 0 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
	weath = 'apr2016'
	mask  = 'new_blu_03'
elif Trim_setting == 'apr16_04':
	id = 1
	wl_lower  = 4700 
	wl_higher = 7200 
	x_srt = 0 	#default 0
	x_end = 920 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
	weath = 'apr2016'
	mask  = 'new_red'
elif Trim_setting == 'apr16_05':
	id = 1
	wl_lower  = 3500 
	wl_higher = 5900 
	x_srt = 0 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
	weath = 'apr2016'
	mask  = 'new_blu'
elif Trim_setting == 'apr16_06':
	id = 1
	wl_lower  = 4700 
	wl_higher = 7200
	x_srt = 0 	#default 0
	x_end = 913 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
	weath = 'apr2016'
	mask  = 'new_red'
elif Trim_setting == 'apr16_07':
	id = 3
	wl_lower  = 4700 
	wl_higher = 7200
	x_srt = 0 	#default 0
	x_end = 912 #default 1024
	y_srt = 0   #default 0
	y_end = 244 #default 245
	weath = 'apr2016'
	mask  = 'new_red_G'
elif Trim_setting == 'mar16_05':
	id = 4
	wl_lower  = 4700 
	wl_higher = 7200
	x_srt = 0   #25 	#default 0
	x_end = 918 #915 #default 1024
	y_srt = 0   #default 0
	y_end = 243 #245 #default 245
	weath = 'mar2016'
elif Trim_setting == 'mar16_06':
	id = 1
	wl_lower  = 3500 
	wl_higher = 5900 
	x_srt = 35 	#default 0
	x_end = 980 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
	weath = 'mar2016'
elif Trim_setting == 'apr15_18':
	id = 1
	wl_lower  = 5000
	wl_higher = 7400
	x_srt = 14 	#default 0
	x_end = 850 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245
	weath = 'apr2015'
	mask  = 'old_red'
elif Trim_setting == 'apr15_16':
	id = 1
	wl_lower  = 3500
	wl_higher = 5800
	x_srt = 35 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 245 #default 245s
	weath = 'apr2015'
	mask  = 'old_blu'
else:
	print 'No trim setting applied'
	id = 1
	wl_lower  = 3500 
	wl_higher = 5900 
	x_srt = 0 	#default 0
	x_end = 1024 #default 1024
	y_srt = 0   #default 0
	y_end = 246 #default 246
	weath = 'apr2016'
	mask  = 'new_blu'

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

filelis = sci_im_list
# for s in sci_im_list:
# 	filename = string.split(s,'/')[-1]
# 	filenum  = filename[-9:-5]
# 	filelis.append(filenum)

print 'FLUX CALIB FOR '+str(name)
print 'Files found: '+str(filelis)

#-------------------------------------#
# FIND WEATHER INFO FROM WEATHER FILE #
#-------------------------------------#

#gives Temp in C, Realtive Humidity, and Pressure in InHg
weather_file = op.join(vp_redux_path, 'weather_files/weather-data_'+weath+'.txt')
print weather_file

#This function is for finding the closest folder date to the data taken
#this is used for finding the best folder to pull special cal data from
def close_datetime(weather_file,data_datetime):
	weather_data = np.genfromtxt(weather_file, dtype=[('date','S10'),('time','S10'),('temp','f8'),('relHum','f8'),('Press','f8')] ,usecols = (1,2,3,4,5), skip_header=1, delimiter=',' , unpack=True)

	weath_datetime_lis = []
	temp   = []
	relHum = []
	Press  = []
	for l in range(len(weather_data)):
		weath_date     = string.split(weather_data[l][0],'-')
		weath_time     = string.split(weather_data[l][1],':')
		weath_datetime = datetime(year=int(weath_date[2]), month=int(weath_date[0]), day=int(weath_date[1]), hour=int(weath_time[0]), minute=int(weath_time[1]))
		weath_datetime_lis.append(weath_datetime)
		temp.append(weather_data[l][2])
		relHum.append(weather_data[l][3])
		Press.append(weather_data[l][4])
	#finds the closest date/folder to the date the data was taken 
	#close_date = min(weath_datetime, key=lambda x:abs(x-data_datetime))
	near_date = min(weath_datetime_lis, key=lambda x: abs(x - data_datetime))
	idx       = weath_datetime_lis.index(near_date)
	clos_temp   = (temp[idx] - 32.0) * (5.0/9.0) #convert Fahrenheit to Celcius
	clos_relHum = relHum[idx]
	clos_Press  = Press[idx]/0.0295300 #convert InHg to mBar (https://www.weather.gov/media/epz/wxcalc/pressureConversion.pdf)

	if np.isnan(clos_temp):
		clos_temp = 7.0
	if np.isnan(clos_relHum):
		clos_relHum = 54.5
	if np.isnan(clos_Press):
		clos_Press = 790.0

	return clos_temp, clos_relHum, clos_Press

#--------------------------------------------#
# FIND OBSERVATION PARAMETERS FROM DITH FILE #
#--------------------------------------------#

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

im1 = pyfits.open(str(filelis[0]))
d1 =  im1[0].data
h1 = im1[0].header

im2 = pyfits.open(str(filelis[1]))
d2 =  im2[0].data
h2 = im2[0].header

im3 = pyfits.open(str(filelis[2]))
d3 =  im3[0].data
h3 = im3[0].header

im4 = pyfits.open(str(filelis[3]))
d4 =  im4[0].data
h4 = im4[0].header

im5 = pyfits.open(str(filelis[4]))
d5 =  im5[0].data
h5 = im5[0].header

im6 = pyfits.open(str(filelis[5]))
d6 =  im6[0].data
h6 = im6[0].header

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

#--------------------------------------------------#
# FIND COORDINATE INFO FROM HEADER AND COORD. FILE #
#--------------------------------------------------#

#Find siderial time from the image header (use this to calculate HA)
h_ST = string.split(h1['ST'],':')
obj_ST = Angle(h_ST[0]+'h'+h_ST[1]+'m'+h_ST[2]+'s')

#Get the objects RA and DEC from the coordiante file 
coord_data = np.genfromtxt(coord_file, dtype=[('obj','S10'),('RA1','S10'),('RA2','S10'),('RA3','S10'),('DEC1','S10'),('DEC2','S10'),('DEC3','S10')] ,usecols = (1,2,3,4,5,6,7), unpack=True)
for o in range(len(coord_data)):
	co = coord_data[o]
	if co[0] == name: 
		obj_RA  = Angle(co[1]+'h'+co[2]+'m'+co[3]+'s')
		obj_DEC = Angle(co[4]+'d'+co[5]+'m'+co[6]+'s')

#calculate the objects HA from ST and RA
obj_HA = obj_ST-obj_RA

print 'Object Coords Found: RA '+str(obj_RA.hms)+' DEC '+str(obj_DEC.dms)+ ' ST: '+str(obj_ST.hms)+' Hour Angle '+str(obj_HA.hms)

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

#---------------------------#
# BUILD WAVELENGTH SOLUTION #
#---------------------------#

wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
start_wave = float(h1['CRVAL1']) #this is the starting wavelength value

#wave_orig = np.add(np.arange(wave_sol,(wave_sol*ln_orig)+(0.5*wave_sol),wave_sol),start_wave)  # wavelength solution for each of the fibers (should be same across the board)
wave_orig = np.add(np.arange(0.0,(wave_sol*(ln_orig-1))+(0.5*wave_sol),wave_sol),start_wave)
wave = wave_orig[x_srt:x_end]	#trim the wavelength solution to make the trimmed image

print 'Data Shape: '+str(np.shape(d1)), 'Wave Sol Shape: '+str(np.shape(wave))

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

	#flux weighted position, output is the positon of the centroid [sx_avg,sy_avg] of star
	sx_avg=sum([wx1,wx2,wx3,wx4,wx5,wx6])/sum([sum(totflux1[sel1]),sum(totflux2[sel2]),sum(totflux3[sel3]),sum(totflux4[sel4]),sum(totflux5[sel5]),sum(totflux6[sel6])])
	sy_avg=sum([wy1,wy2,wy3,wy4,wy5,wy6])/sum([sum(totflux1[sel1]),sum(totflux2[sel2]),sum(totflux3[sel3]),sum(totflux4[sel4]),sum(totflux5[sel5]),sum(totflux6[sel6])])

	if np.isnan(sx_avg) or np.isnan(sy_avg):
		sys.exit('CAN\'T FIND THE STAR: please check the central fiber defined for your standard star (fib0)')

	print 'position of centroid of star: '
	print '	x: '+str(sx_avg)
	print '	y: '+str(sy_avg)

 	#--------------------------------------#
	# PLOT THE SUMMED SPECTRUM OF THE STAR #
	#--------------------------------------#
	totflux1 = sum(aux1, axis = 0)
	totflux2 = sum(aux2, axis = 0)
	totflux3 = sum(aux3, axis = 0)
	totflux4 = sum(aux4, axis = 0)
	totflux5 = sum(aux5, axis = 0)
	totflux6 = sum(aux6, axis = 0)

	plt.plot(wave,totflux1,label='d1')
	plt.plot(wave,totflux2,label='d2')
	plt.plot(wave,totflux3,label='d3')
	plt.plot(wave,totflux4,label='d4')
	plt.plot(wave,totflux5,label='d5')
	plt.plot(wave,totflux6,label='d6')
	plt.ylabel('Relative Flux')
	plt.xlabel('Wavelength (A)')
	plt.title('Total Summed Spec')
	plt.legend()
	plt.savefig('Summed_Standard_Spec_'+str(name)+'.png')
	#plt.show()

	#-------------------------------------------------------#
	# ATMOSPHERIC DIFFERENTIAL REFRACTION CORRECTION - CURE #
	#-------------------------------------------------------#

	#find weather information by finding time and date from header of 3rd dither
	#example header format: 2017-03-22T18:33:26
	data_date = string.split(h3['date'][0:10],'-')
	data_time = string.split(h3['date'][11:20],':')
	data_datetime = datetime(year=int(data_date[0]), month=int(data_date[1]), day=int(data_date[2]), hour=int(data_time[0]), minute=int(data_time[1]), second=int(data_time[2])) 
	clos_temp, clos_relHum, clos_Press = close_datetime(weather_file,data_datetime)

	print 'Weather: temp:'+str(clos_temp), 'RH: '+str(clos_relHum), 'Press: '+str(clos_Press)

	#takes in astropy Angle objects, returns angle in degrees
	def slit_offset(Ha,lat,dec):
		Ha  = Ha.radian
		dec = dec.radian
		lat = math.radians(lat)

		num = math.sin(Ha)*math.cos(lat)
		mid_den = ((math.sin(lat)*math.sin(dec)) + (math.cos(lat)*math.cos(dec)*math.cos(Ha)))**2
		den = (1-mid_den)**(0.5)
		off_angle = math.asin(num/den)
		return math.degrees(off_angle)

	# float Projector::atmos_refraction( float l, float secz ) const
	# Temperature [C] # Relative Humidity [%] # Pressure [mbar]
	# returns arcsecond offset 
	def atmos_refraction(secz,ref_wl, wl, TC = 7.0, RH = 54.5, P = 600.0):
	#def atmos_refraction(wl, coef, secz, ref_wl, TC = 7.0, RH = 54.5, P = 600.0):

		ref_wl = ref_wl/1e4
		wl = np.divide(wl,1e4)

		T=TC+273.16

		PS=-10474.0+116.43*T-0.43284*T*T+0.00053840*T*T*T

		P2=RH/100.0*PS
		P1=P-P2

		D1=P1/T*(1.0+P1*(57.90*1.0E-8-(9.3250*1.0E-4/T)+(0.25844/(T*T))));
		D2=P2/T*(1.0+P2*(1.0+3.7E-4*P2)*(-2.37321E-3+(2.23366/(T*T))-(710.792/(T*T))+(7.75141E4/(T*T*T))))

		N0_1 = nl(ref_wl, D1, D2)
		N1_1 = nl(wl, D1, D2)

		return tan(arccos(1.0/secz))*(N0_1-N1_1)*206264.8

	#float Projector::nl( float l, float D1, float D2 ) const
	def nl(wl, D1, D2):
		S0 = np.divide(1,wl)
		S2 = np.multiply(S0,S0)
		S4 = np.multiply(S2,S2)
		S6 = np.multiply(S4,S2)
		#N_1 = 1.0E-8*((2371.34+683939.7/(130.0-S2)+4547.3/(38.9-S2))*D1+(6487.31+58.058*S2-0.71150*S4+0.08851*S6)*D2)
		N_1 = 1.0E-8*((2371.34+480000.0/(130.0-S2)+4547.3/(38.9-S2))*D1+(6487.31+58.058*S2-0.71150*S4+0.08851*S6)*D2) #coef=480000.0
		return N_1

	#shift x and y positions for 
	#only shifts in x direction assuming the field is lined up perp. to horizon. 
	#ADR_shift = atmos_refraction(secz= airmass ,ref_wl= wave[0], wl= wave, TC = 7.0, RH = 54.5, P = 600.0)
	print 'Reference Wavelength: '+str(wave[(len(wave)/2)]), str(len(wave)/2)+'/'+str(len(wave))
	ADR_mag   = atmos_refraction(secz= airmass ,ref_wl= wave[(len(wave)/2)], wl= wave, TC = clos_temp, RH = clos_relHum, P = clos_Press)

	off_angle = slit_offset(obj_HA, McD_lat, obj_DEC)
	print 'Parallactic Angle: '+str(off_angle)

	#------------------------------------------------------#
	# ATMOSPHERIC DIFFERENTIAL REFRACTION CORRECTION - COM #
	#------------------------------------------------------#

	print '[7] FINDING CENTOID OFFSETS FOR ADR CORRECTION'

	#this will find the ADR correction by finding calculating the center of the star at every wavelength to calculate the shift in the star due to ADR effects
	#to attempt to somwhat 'smooth' the solution it will caluclate the center with a few wavelenths (pixels) on either side of the one being used. 

	#this is the factor that you want to smooth it by. 
	#will include +/- sfac pixels when calculating the center for that pixel (wavelength)
	sfac = 50
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

		#center of mass y value
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
	ADR_x_norm = np.subtract(ADR_x, ADR_x[(len(wave)/2)])
	ADR_y_norm = np.subtract(ADR_y, ADR_y[(len(wave)/2)])

	def line(x, m, b):
		return (m*x)+b

	def angle_from_m(x, m, b):
		y1 = line(x[0], m, b)
		y2 = line(x[-1], m, b)
		m_to_ang = math.atan((y2-y1)/(x[-1]-x[0]))
		return math.degrees(m_to_ang)

	#need to give it the angle of the cos (90-para_angle)
	def corr_mag_ang(ang, ADR_mag):
		cos_off =  math.cos(math.radians(ang))
		ADR_shift = ADR_mag/cos_off
		return ADR_shift

	def rotate_line(origin, x, y, angle):
	    """
	    Rotate a point counterclockwise by a given angle around a given origin.
	    """
	    ox, oy = origin

	    #converts user input degrees to radians
	    angle = np.deg2rad(angle)

	    qx = ox + np.multiply(math.cos(angle), np.subtract(x, ox)) - np.multiply(math.sin(angle), np.subtract(y, oy))
	    qy = oy + np.multiply(math.sin(angle), np.subtract(x, ox)) + np.multiply(math.cos(angle), np.subtract(y, oy))

	    return qx, qy

	#fit line to the data and find the offset angle 
	popt, pcov = curve_fit(line, ADR_x, ADR_y) #, sigma = sigma, absolute_sigma=True)
	m_ang = angle_from_m(ADR_x, popt[0], popt[1])
	para_angle = 90-m_ang
	print 'Parallactic Angle from fit: '+str(para_angle)

	#find the shift from the reference point to the fit of the reference point
	y_orgin_shift = ADR_y[(len(wave)/2)] - line(ADR_x[(len(wave)/2)],popt[0],popt[1])
	origin = (ADR_x[(len(wave)/2)], ADR_y[(len(wave)/2)] - y_orgin_shift)
	print 'Origin: '+str(origin)

	#Correct for the magnitiude of the line based on the angle found from data
	ADR_shift = corr_mag_ang(m_ang,ADR_mag)
	ADR_shift = np.add(ADR_shift, origin[1])
	fit_ADR_x, fit_ADR_y = rotate_line(origin, np.add(np.zeros(len(ADR_shift)),origin[0]), ADR_shift , -(180+para_angle))

	#Shift assuming the calculated parallactic angle 
	shift_xc, shift_yc = rotate_line(origin, np.add(np.zeros(len(ADR_shift)),origin[0]), ADR_shift , -(180+off_angle))

	poly_x = np.polyfit(wave, ADR_x, 2)
	poly_y = np.polyfit(wave, ADR_y, 2)

	px = np.poly1d(poly_x)
	py = np.poly1d(poly_y)

	fit_x = px(wave)
	fit_y = py(wave)

	print 'Airmass: '+str(airmass)

	#plot the shift in the centriod as a function of wavelength 
	fig = plt.figure(figsize=(7,9))
	ax1 = fig.add_subplot(311)

	#plot circle at the found center of the star 
	circle1=plt.Circle(origin, 0.002, color='grey',label='origin')
	fig.gca().add_artist(circle1)

	#plot the shift in the centriod as a function of wavelength 
	MAP ='rainbow_r'
	cm  = plt.get_cmap(MAP)

	#plot the shift from data found from COM calculation
	ax1.scatter(ADR_x, ADR_y, s=10, c = range(len(ADR_x)), cmap = cm, label='Data COM')

	#plot the ADR_shift reference point
	ax1.scatter(np.add(np.zeros(len(ADR_shift)),origin[0]), ADR_shift, s=2, color='black', label='shift ref')
	#plot the line fit to that data 
	ax1.scatter(ADR_x,line(ADR_x, popt[0], popt[1]), s=2, color='blue', label='fit')
	#plot the ADR_shifted line from fitted parallactic angle
	ax1.scatter(fit_ADR_x, fit_ADR_y, s=2, color='green', label='fit PA')
	#plot the ADR_shifted line from calculated parallactic angle 
	ax1.scatter(shift_xc,shift_yc, s=2, color='red', label='calc PA')

	ax1.scatter(fit_x,fit_y, s=2, color='grey',label='fit')

	ax1.set_xlabel('X-Shift (arcsec)')
	ax1.set_ylabel('Y-Shift (arcsec)')
	plt.legend()

	#plt.text(0.8,0.45,'C_PA: '+str(round(para_angle,2)),transform = ax1.transAxes, color='navy',size='small', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))		
	#plt.text(0.8,0.35,'H_PA: '+str(round(off_angle,2)),transform = ax1.transAxes, color='navy',size='small', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))		
	#plt.text(0.8,0.25,'Temp: '+str(round(clos_temp,2)),transform = ax1.transAxes, color='navy',size='small', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))		
	#plt.text(0.8,0.15,'Press: '+str(round(clos_Press,2)),transform = ax1.transAxes, color='navy',size='small', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))		
	#plt.text(0.8,0.05,'relHum: '+str(round(clos_relHum,2)),transform = ax1.transAxes, color='navy',size='small', bbox=dict(facecolor='none', edgecolor='navy', pad=10.0))		

	ax2 = fig.add_subplot(312)
	ax2.scatter(wave, fit_ADR_x, s=1, color='red', label='x_fit')
	ax2.scatter(wave, ADR_x, s=1, color='orange', label='x_data')
	ax2.scatter(wave, fit_x, s=1, color='green', label='fit')
	ax2.set_xlabel('wavelength (A)')
	ax2.set_ylabel('X-Shift (arcsec)')
	plt.legend()

	ax3 = fig.add_subplot(313)
	ax3.scatter(wave, fit_ADR_y, s=1, color='blue', label='y_fit')
	ax3.scatter(wave, ADR_y, s=1, color='purple', label='y_data')
	ax3.scatter(wave, fit_y, s=1, color='green', label='fit')
	ax3.set_xlabel('wavelength (A)')
	ax3.set_ylabel('Y-Shift (arcsec)')
	plt.legend()

	fig.savefig('ADR_StarCenShift_'+str(name)+'.png')
	plt.show()

	#------------------------------------#
	# MAKE 2D GAUSSIAN MODEL OF THE STAR #
	#------------------------------------#

	print '[8] BUILDING MODELS OF STANDARD STAR'

	Npts    = 1e5 # should be enough points
	tot_fibs = y_end*len(filelis) #total number of fibers = number of fibers times number of dithers

	fiber_x = np.hstack((fx1,fx1,fx2,fx4,fx5,fx6)) # list of fiber x postions from all dithers
	fiber_y = np.hstack((fy1,fy1,fy2,fy4,fy5,fy6)) # list of fiber y postions from all dithers

	fiber_flux   = np.vstack((aux1,aux2,aux3,aux4,aux5,aux6)) # This is where all of the fiber fluxes are kept ()
	F            = np.zeros((len(wave),)) # least squares method

	for j in range(len(wave)):

		if j%10 == 0:
			print 'Model '+str(j)+' of '+str(len(wave))

		mean = (fit_x[j],fit_y[j])
		cov = seeing * np.eye(2)
		P = np.random.multivariate_normal(mean,cov, int(Npts))
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

	#stand_spec = stand_spec[x_srt:x_end]

	#get the wavelength values for this spectrum 
	wave_sol = float(h1['CDELT1']) #this is the wavelength per pixel
	start_wave = float(h1['CRVAL1']) #this is the starting wavelength value
	wave_orig = np.add(np.arange(0.0,(wave_sol*(ln_orig-1))+(0.5*wave_sol),wave_sol),start_wave)
	wave = wave_orig[x_srt:x_end] 	#trim the wavelength solution to make the trimmed image

	#np.divide by exposure time to get per second
	stand_spec_perpix = np.divide(stand_spec,exptime) #DN/s/pixel
	#convert standard spectrum to be per anstrom so it is invariant
	stand_spec = np.divide(stand_spec_perpix,wave_sol) #DN/s/A

	#--------------------------------#
	# INTERPOLATE OVER KNOW SPECTRUM #
	#--------------------------------#

	#interpolate the know spectrum and create the spectrum over the wavelength array for the standard star spectrum
	f = interpolate.interp1d(x=sl, y=sflam)
	sflamint = f(wave)

	#-------------------------------#
	# DEFINE SPECTRAL FEATURES MASK #
	#-------------------------------#	

	def near_idx(a,val):
	    idx = (np.abs(a-val)).argmin()
	    return idx

	#define regions of the spectrum that are featureless 
	if mask == 'old_red':
		feature_mask = [[wave[0],5280], [5340,5400], [5430,5550], [5590,5830], [5930,6275], [6380,6530], [6680,6860], [6890,wave[-1]]]
	if mask == 'old_blu':
		feature_mask = [[wave[0],3880], [3910,3960], [3990,4080], [4120,4190], [4240,4330], [4380,4670], [4700,4825], [4880,5400], [5440,5560], [5600,wave[-1]]]
	if mask == 'new_red':
		feature_mask = [[wave[0],4840], [4900,5280], [5350,5400], [5420,6270],[6290,6550], [6590,6870], [6890,wave[-1]]]
	if mask == 'new_red_G':
		feature_mask = [[wave[0],4830], [4950,5280], [5350,5400], [5420,5850], [5950,6270],[6290,6550], [6660,6870], [6890,wave[-1]]]
	if mask == 'new_blu_03':
		feature_mask = [[wave[0],3830], [3850,3880], [3910,3960], [3990,4080], [4120,4190], [4240,4250], [4290,4330], [4380,4670], [4700,4825], [4880,4940], [4980,5320], [5360, 5400], [5440,wave[-1]]]
	if mask == 'new_blu':
		feature_mask = [[wave[0],3830], [3850,3880], [3910,3960], [3990,4080], [4120,4190], [4240,4250], [4290,4330], [4380,4400], [4600,4670], [4700,4825], [4880,5320], [5360, 5400], [5440,wave[-1]]]


	stand_spec_sm = []
	known_spec_sm = []
	wave_sm       = []

	for i in range(len(feature_mask)):
		sect = feature_mask[i]
		up_idx = near_idx(wave,sect[1])
		lo_idx = near_idx(wave,sect[0])
		sm_scale = 9
		mov_idx = lo_idx

		for j in range(((up_idx - lo_idx)/(sm_scale + 1)) - 1):
			avg_wave   = (wave[mov_idx+sm_scale] + wave[mov_idx])/2.0
			st_spec_pt = np.median(stand_spec[mov_idx:mov_idx+sm_scale])
			kn_spec_pt = np.median(sflamint[mov_idx:mov_idx+sm_scale])
			
			mov_idx = mov_idx + sm_scale + 1

			stand_spec_sm.append(st_spec_pt)
			known_spec_sm.append(kn_spec_pt)
			wave_sm.append(avg_wave)

	#-------------------------------#
	# CALCULATE SENISIVITY FUNCTION #
	#-------------------------------#

	print '[9] FINDING SENSITIVITY FUNCTION'

	#np.divide to get the sensitivity funciton 
	sens_func = np.divide(stand_spec, sflamint) #DN*cm^2 / erg
	sens_func_sm = np.divide(stand_spec_sm, known_spec_sm) #smoothed version 
	spl = UnivariateSpline(wave_sm, sens_func_sm,ext=3)
	sens_func_spline = spl(wave)
	
	np.save('SensitivityFunc_'+str(name), sens_func_spline)
	#np.save('WavelengthSolu_'+str(name), wave)

	#print sens_func_spline

	#------------#
	# PLOT STUFF #
	#------------#

	plt.subplot(4,1,1)
	plt.plot(wave, stand_spec, color='blue', label='Standard')
	plt.plot(wave_sm, stand_spec_sm, color = 'red', ls='None', marker = '.', ms = 7, label='Stand Sm')
	plt.title('Observed Spectrum')
	plt.ylabel('DN/s/pixel')
	plt.xlabel('Wavelength (A)')
	#plt.legend()

	plt.subplot(4,1,2)
	plt.plot(wave, sflamint, color = 'red', label='Known Spec')
	plt.plot(wave_sm, known_spec_sm, color = 'purple', ls='None', marker = '.', ms = 7, label='Known Sm')
	plt.title('Oke Spectrum')
	plt.ylabel('Flux (ergs/s/cm^2/A)')
	plt.xlabel('Wavelength (A)')
	#plt.legend()

	plt.subplot(4,1,3)
	plt.plot(wave, sens_func, color = 'purple', label='Known Sm')
	plt.plot(wave_sm, sens_func_sm, color = 'orange', ls='None', marker = '.', ms = 7, label='Known Sm')
	plt.plot(wave, sens_func_spline, color = 'black', label='Known Sm')
	plt.title('Sensitivity Function (Observed/Oke)')
	plt.ylabel('DN*cm^2 / erg')
	plt.xlabel('Wavelength (A)')
	#plt.legend()
	#plt.show()

	plt.subplot(4,1,4)
	corr_stand_spec = stand_spec/sens_func_spline
	plt.plot(wave, corr_stand_spec, color = 'green', label='Known Sm')
	plt.title('Standard Spectrum divided by Sensitivity Function')
	plt.ylabel('Ergs/sec/cm^2/A')
	plt.xlabel('Wavelength (A)')

	#plt.tight_layout()

	plt.savefig('Flux_Calib_'+str(name)+'_'+Trim_setting+'.png')
	plt.show()



