import string
import subprocess
import os.path as op
import os 
import sys

#set up path to your cure bin for NEW CURE
curebin_path  = '/Users/Briana/Documents/cure_v2/cure/bin'
cen_file = op.join(curebin_path,'../config/IFUcen_VP2_27m.txt')

###################
# Setting CUREBIN #
###################

#checking that LRS2 is defined in specconf.h 
cureversion = os.popen(op.join(curebin_path, 'cureversion')).readlines()
spec_define = cureversion[4].split(' ')[1]
instrument = spec_define.rstrip('\n')

if instrument == 'VIRUS-P':
    print ('NEW CURE is set for VIRUS-P reduction')
else: 
    print ('You need to update specconf.h in CURE to define VIRUS-P')
    sys.exit('Right now specconf.h defines '+instrument)

#make cube parameters

#Regridded sample size on sky in arcsec 
skysamp = 2.0 #chosen to be 2 becuase have the fiber size to properly nyquist sample 
#Samples further away will haver weight=0 [arcsec]
maxdist = 15  #chose to be a little over 3x the size of the fiber
#Gaussian sigma of interpolation kernel [arcsec]. 
sigma   = 1.8 #This was chosen using FWHM = 2*sqt(2ln(2))*sigma. use FWHM = 4.24 (size of fiber). sigma = (4.24/2.355)

#find the filename numbers from sci_im_lis.lis
sci_im_list = open('sci_im_lis.lis').read().splitlines()
filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

print 'Files found: '+str(filelis)

for i in range(len(filelis)):
	bash_com = op.join(curebin_path,'mkcube -r Fcal -i '+cen_file+' -a '+str(skysamp)+' -k '+str(maxdist)+' -s '+str(sigma)+' dith2.txt')
	print bash_com
	subprocess.call(bash_com, shell=True)