import string
import subprocess
import os.path as op
import os 
import sys
from astropy.io import fits
import numpy as np

#run this inside of the object folder you want to build a data cube of

#************************#
# User Defined Variables #
#************************#

#set up path to your cure bin for NEW CURE
curebin_path  = '/Users/Briana/Documents/cure_v2/cure/bin'
vp_redux_path = '/Users/Briana/Documents/Grad_School/VIRUS_P/VP_reduction'

#If your headers are missing RA and DEC can read in your work file to add them to headers
add_coords = True
#you only need this path if you want to add RA and DEC to headers
worklist = '/Users/Briana/Documents/Grad_School/VIRUS_P/Observing/mar2017_run/object_lis.txt'

#----------------------#
# make cube parameters #
#----------------------#

#Regridded sample size on sky in arcsec 
skysamp = 2.0 #chosen to be 2 becuase have the fiber size to properly nyquist sample 
#Samples further away will haver weight=0 [arcsec]
maxdist = 15  #chose to be a little over 3x the size of the fiber
#Gaussian sigma of interpolation kernel [arcsec]. 
sigma   = 1.8 #This was chosen using FWHM = 2*sqt(2ln(2))*sigma. use FWHM = 4.24 (size of fiber). sigma = (4.24/2.355)

#+++++++++++++++++++ end user defined variables +++++++++++++++++++++#

#*****************#
# Setting CUREBIN #
#*****************#

#checking that LRS2 is defined in specconf.h 
cureversion = os.popen(op.join(curebin_path, 'cureversion')).readlines()
spec_define = cureversion[4].split(' ')[1]
instrument = spec_define.rstrip('\n')

if instrument == 'VIRUS-P':
    print ('NEW CURE is set for VIRUS-P reduction')
else: 
    print ('You need to update specconf.h in CURE to define VIRUS-P')
    sys.exit('Right now specconf.h defines '+instrument)

cen_file = op.join(vp_redux_path,'IFUcen_VP2_27m_mkcube.txt')

#find the filename numbers from sci_im_lis.lis
sci_im_list = open('sci_im_lis.lis').read().splitlines()

#***************************************#
# If add_coords: Reading in Coordinates #
#***************************************#

if add_coords:
	#read in the worklist 
	work_data = np.genfromtxt(worklist, usecols = (1,2,3,4,5,6,7), dtype=None)
	#open the first image to get header information 
	hdr = fits.getheader('FcalFe'+string.split(sci_im_list[0],'/')[-1], 0)
	#This is the object name in the header without dither num
	obj = string.split(hdr['OBJECT'],' ')[0]

	#find the line in work list that has the same object name and save the coordinates of that object 
	work_obj = None
	for d in work_data:
		if d[0] == obj:
			work_obj = d
	if work_obj == None:
		sys.exit('Your object was not found in the work list')
	else:
		print 'Coordinates added will be added to FeCal headers:'
		print '	RA', str(work_obj[1])+':'+str(work_obj[2])+':'+str(work_obj[3])
		print '	DEC', str(work_obj[4])+':'+str(work_obj[5])+':'+str(work_obj[6])

#***********************************************#
# Finding Files and Adding Coords if add_coords #
#***********************************************#

filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

	#adding the coordinates to the headers id add_coords
	if add_coords:
		dat, hdr = fits.getdata('FcalFe'+filename, header=True)
		hdr.set('RA      ', str(work_obj[1])+':'+str(work_obj[2])+':'+str(work_obj[3]))
		hdr.set('DEC      ', str(work_obj[4])+':'+str(work_obj[5])+':'+str(work_obj[6]))
		fits.writeto('FcalFe'+filename, dat, hdr, clobber=True)

#Checks for missing dithers before bulding data cube 
if len(filelis) == 3 or 6:
	print 'Found '+str(len(filelis))+' Dithers: '+str(filelis)
else: 
	print 'ONLY Found '+str(len(filelis))+' Dithers: '+str(filelis)
	sys.exit('Check for Missing Files')

#####################
# Run CURE's mkcube #
#####################

bash_com = op.join(curebin_path,'mkcube -r FcalFe -i '+cen_file+' -a '+str(skysamp)+' -k '+str(maxdist)+' -s '+str(sigma)+' dith2.txt')
print bash_com
subprocess.call(bash_com, shell=True)

