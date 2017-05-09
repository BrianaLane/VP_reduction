import string
import subprocess
import glob
import os.path as op
import shutil 

#************************#
# User Defined Variables #
#************************#

#set up path to your cure bin
curebin_path  = '/Users/Briana/Documents/cure/virusp1/bin'

#prefix for your sciene files 
prefix = 'vp'

#extract error files
extract_error = True

#define wavelength range (in angstroms)
wl_lower = 3500
wl_higher = 5900

#+++++++++++++++++++ end user defined variables +++++++++++++++++++++#

#*********************************#
# Find files and run fiberextract #
#*********************************#

sci_im_list = glob.glob('Sp*.fits')	

filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

print 'Files found: '+str(filelis)
for i in range(len(filelis)):
	if extract_error:
		bash_com = op.join(curebin_path,'fiberextract -d pes'+str(prefix)+str(filelis[i])+'.dist -f pes'+str(prefix)+str(filelis[i])+'.fmod -c -r 1,246 -l '+str(wl_lower)+','+str(wl_higher)+' e.Spes'+str(prefix)+str(filelis[i])+'.fits')
	else:
		bash_com = op.join(curebin_path,'fiberextract -d pes'+str(prefix)+str(filelis[i])+'.dist -f pes'+str(prefix)+str(filelis[i])+'.fmod -c -r 1,246 -l '+str(wl_lower)+','+str(wl_higher)+' Spes'+str(prefix)+str(filelis[i])+'.fits')
	#$CUREBIN/fiberextract -d pesvp0068.dist -f pesvp0068.fmod -c -r 1,246 -l [3500,5800] Spesvp0068.fits 
	subprocess.call(bash_com, shell=True)

if extract_error:
	for i in range(len(filelis)):
		shutil.move('Fee.Spes'+str(prefix)+str(filelis[i])+'.fits','e.FeSpes'+str(prefix)+str(filelis[i])+'.fits')
