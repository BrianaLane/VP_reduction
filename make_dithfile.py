import numpy as np
import pyfits 
import string
import glob
import subprocess
import shutil
import os.path as op
import os
import warnings
import sys

#This script should be run inside the folder of the science object you are working on
#You should have already run match_guiderFrames.py and this folder should contain sci_im_list.lis
#This script calls source extractor so you must have it installed
#	You must set the param file to print out flux and FWHM_IMAGE

#************************#
# User Defined Variables #
#************************#

#define the field youa are working on. 
#the script will use the guide start coordinates defined for that field
field = 'F4'

#This path should point to the location of both files:
#	runsex.in, dith_vp2_subdither.lis
script_path = '/Users/Briana/Documents/Grad_School/VIRUS_P/VP_reduction'

#the number (starting from zero) of the column for flux in .cat files
flux_num = 1
#the number (starting from zero) of the column for the fwhm (pixels) in the .cat files
fwhm_num = 7

#***********************************************#
# User defined guide star corrdinates and field #
#***********************************************#

#below you should define the corrdinates of the guide star for each field

if field == 'Feige34':
	x1 = 260 #min
	x2 = 300 #max
	y1 = 380 #min
	y2 = 410 #max

elif field == 'GRW70d5824':
	x1 = 360 #min
	x2 = 400 #max
	y1 = 325 #min
	y2 = 355 #max

elif field == 'BD_75D325':
	x1 = 365 #min
	x2 = 405 #max
	y1 = 265 #min
	y2 = 305 #max

elif field == 'Fa':
	x1 = 295 #min
	x2 = 335 #max
	y1 = 380 #min
	y2 = 420 #max

elif field == 'Fb': 
	x1 = 410 #min
	x2 = 470 #max
	y1 = 190 #min
	y2 = 250 #max

elif field == 'Fc':
	x1 = 20 #min
	x2 = 65 #max
	y1 = 230 #min
	y2 = 280 #max

elif field == 'Fi':
	x1 = 100 #min
	x2 = 150 #max
	y1 = 380 #min
	y2 = 430 #max

elif field == 'F4':
	x1 = 100 #min
	x2 = 170 #max
	y1 = 265 #min
	y2 = 330 #max

elif field == 'F7':
	x1 = 330 #min
	x2 = 355 #max
	y1 = 225 #min
	y2 = 250 #max

elif field == 'F8':
	x1 = 80 #min
	x2 = 130 #max
	y1 = 300 #min
	y2 = 340 #max

elif field == 'F9':
	x1 = 300 #min
	x2 = 340 #max
	y1 = 130 #min
	y2 = 170 #max

elif field == 'F11':
	x1 = 275 #min
	x2 = 300 #max
	y1 = 145 #min
	y2 = 170 #max

elif field == 'F10':
	x1 = 110 #min
	x2 = 150 #max
	y1 = 110 #min
	y2 = 145 #max

elif field == 'F12':
	x1 = 465 #min
	x2 = 500 #max
	y1 = 420 #min
	y2 = 465 #max

elif field == 'F13':
	x1 = 280 #min
	x2 = 310 #max
	y1 = 365 #min
	y2 = 400 #max

elif field == 'F14':
	x1 = 430 #min
	x2 = 470 #max
	y1 = 230 #min
	y2 = 270 #max

elif field == 'F15':
	x1 = 170 #min
	x2 = 250 #max
	y1 = 95 #min
	y2 = 170 #max

#********************************#
# Building the guide star frames #
#********************************#

else:
	sys.exit("You did not pick a field with defined corrdinates")

#reading the list of science images
sci_im_list = open('sci_im_lis.lis').read().splitlines()

#make directory to store guider_small frames 
if os.path.exists('guider_small'):
    shutil.rmtree('guider_small')
os.makedirs('guider_small')

#Iterating through science frame list becuase that is how many time_lists there are
for s in range(len(sci_im_list)):
	print 'Science Frame: '+sci_im_list[s]

	#opening the time list for this frame
	dat_name = 'time_'+str(s)+'.lis'
	print 'Time List: '+dat_name
	dat_lis = open(str(dat_name)).read().splitlines()

	#iterating though each guider frame in time list
	for d in dat_lis:
		print "	"+d

		#defining the name of the output file
		origname = string.split(d,'/')[-1] 
		outname  = 'guider_small/'+string.split(origname,'.')[0]+'small.fits'

		#opening the guider frame file
		im = pyfits.open(d)
		dat =  im[0].data

    	#defining small array around guide star in guider frame based on user defined coordiantes 
		small_im = dat[y1:y2,x1:x2]
    	#writing new fits file containing just small frame of guide star
		hdu = pyfits.PrimaryHDU(small_im)
		hdu.writeto(outname, clobber=True)

	print '\n'

#*******************************************#
# Run Source Extractor on small star frames #
#*******************************************#

#copy runsex.in to your directory
#this file is necessary for sorce extractory
shutil.copy ( op.join(script_path,'runsex.in'), './' )

#These commands set up the small.fits guider files for source extractor
#they then run sorce extractor on the files
os.system('ls guider_small/guider*small.fits > small.list')
os.system('awk \'{print "sex",$1,"-c runsex.in -CATALOG_NAME",$1 "cat"}\' small.list > small.tmp')
os.system('sed \'s/fitscat/cat/g\' small.tmp > sexcommand')
os.system('chmod u=rwx sexcommand')
os.system('./sexcommand')

#****************************************************#
# Find Values for Seeing, Airmass, and Relative Flux #
#****************************************************#

flux_avg_lis = []
airm_avg_lis = []
fwhm_avg_lis = []

#iterate through each scienc image 
#calculates seeing,airmass, and relative flux for each image
for s in range(len(sci_im_list)):
	#open time list of guider files for that frame
	dat_name = 'time_'+str(s)+'.lis'
	print "Finding values for guider files in: "+dat_name
	dat_lis = open(dat_name).read().splitlines()

	flux_lis = []
	fwhm_lis = []
	airm_lis = []
	#for each image it iterates through the guider frames to get values for seeing, airmass,and relative flux
	for d in range(len(dat_lis)):
		#for each guider frame open sorce extractor info file 
		origname = string.split(dat_lis[d],'/')[-1] 
		namein  = 'guider_small/'+string.split(origname,'.')[0]+'small.cat'
		print namein
		#from the file read the flux and fwhm (seeing) of the star
		#the try except catches if the file is empty (did not find a star in that frame)
		try:
			flux_aper, fwhm_pix = np.loadtxt(namein, usecols = (flux_num,fwhm_num), unpack=True)
			#this if statement checks if it found more than one star in the file based on if the output is an array
			#if it finds more than one star in a file it prints the output but does not use this file
			if type(flux_aper) is np.ndarray:
				print namein+" file has more than one star"
			else:
				flux_lis.append(flux_aper) # this later gets normalized to the units don't matter
				fwhm = fwhm_pix * 0.51 # 0.51 arcseconds per pixel
				fwhm_lis.append(fwhm)
		except ValueError:
			print namein+" file is empty."

		#open the guider frame to get the airmass from the header
		im = pyfits.open(dat_lis[d])
		hdr =  im[0].header
		airmass = hdr['AIRMASS']
		airm_lis.append(float(airmass))

	#for each science frame append the average of the valuse from all guider frames for that image
	fwhm_avg_lis.append(np.average(fwhm_lis))
	airm_avg_lis.append(np.average(airm_lis))
	flux_avg_lis.append(np.average(flux_lis))

	print '\n'

#normalize all of the flux values the science frames by the max flux value 
max_flux = np.max(flux_avg_lis)
normflux_avg_lis = np.divide(flux_avg_lis,max_flux)

#*************************************#
# Write Dither File With These Values #
#*************************************#

#load the dither list that has the dither pattern for vp2_subdither
x, y = np.loadtxt(op.join(script_path,'dith_vp2_subdither.lis'), usecols = (0,1), unpack=True)

f1 = open('dith.txt','w')
f2 = open('dith2.txt','w')

for s in range(len(sci_im_list)):
	sci = string.split(sci_im_list[s],'/')[-1]
	pre = string.split(sci,'S')[1]
	prefix = string.split(pre,'.')[0]
	print 'Adding '+prefix+' to dith file'

	files = 'S'+str(prefix)+'.fits'
	pmod = str(prefix)+'.pmod'
	dist = str(prefix)+'.dist'
	fmod = str(prefix)+'.fmod'

	x1 = x[s]
	y1 = y[s]
	see  = float(fwhm_avg_lis[s])
	norm = float(normflux_avg_lis[s])
	air  = float(airm_avg_lis[s])

	f1.write(str(files)+'	'+str(pmod)+'	'+str(dist)+'	'+str(fmod)+'	'+str(x1)+'	'+str(y1)+'	'+str(see)+'	'+str(norm)+'	'+str(air)+'\n') 
	f2.write(str('S'+str(prefix))+'	'+str(prefix)+'	'+str(x1)+'	'+str(y1)+'	'+str(see)+'	'+str(norm)+'	'+str(air)+'\n') 

f1.close()
f2.close()

print('your dither.txt file has been created')
print('BOOM')