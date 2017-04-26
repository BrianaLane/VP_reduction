import string
import subprocess
import os.path as op

#set up path to your cure bin
curebin_path  = '/Users/Briana/Documents/cure/virusp1/bin'

#find the filename numbers from sci_im_lis.lis
sci_im_list = open('sci_im_lis.lis').read().splitlines()
filelis = []
for s in sci_im_list:
	filename = string.split(s,'/')[-1]
	filenum  = filename[-9:-5]
	filelis.append(filenum)

print 'Files found: '+str(filelis)

fiberextract=True

#define wavelength range (in angstroms)
wl_lower = 3500
wl_higher = 5900

fileprefix = 'Spesvp'
DMprefix = 'pesvp'

if fiberextract:
	print '[3] PERFORMING FIBER EXTRACTION'
	for i in range(len(filelis)):
		bash_com = op.join(curebin_path,'fiberextract -d '+str(DMprefix)+str(filelis[i])+'.dist -f '+str(DMprefix)+str(filelis[i])+'.fmod -c -r 1,246 -l '+str(wl_lower)+','+str(wl_higher)+' '+str(fileprefix)+str(filelis[i])+'.fits')
		#$CUREBIN/fiberextract -d pesvp0068.dist -f pesvp0068.fmod -c -r 1,246 -l [3500,5800] Spesvp0068.fits 
		subprocess.call(bash_com, shell=True)