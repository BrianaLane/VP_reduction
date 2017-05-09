import numpy as np
import pyfits
import glob
import os.path as op
import string
import cosmics

standard_stars = ['Feige34','Feige67','GRW70d5824','BD_75d325','BD_75D325']

#this function calls cosmics to do cosmic ray rejection on a list of filenames 
def rmcosmicfits(filenames):
    
    for i in xrange(len(filenames)):

        print (filenames[i])

        array, header = cosmics.fromfits(filenames[i])
        # array is a 2D numpy array

        im_gain = header['GAIN1']
        im_RN = header['RDNOISE1']

        # Build the object :
        c = cosmics.cosmicsimage(array, gain=im_gain, readnoise=im_RN, satlevel = 65535.0, sigclip = 5.0, sigfrac = 0.3, objlim = 7.0)
        # There are other options, check the manual...

        # Run the full artillery :
        c.run(maxiter = 4)

        # Write the cleaned image into a new FITS file, conserving the original header :
        cosmics.tofits(filenames[i], c.cleanarray, header)

        # If you want the mask, here it is :
        #cosmics.tofits("s20160909T093737.7_066LL_sci_mask.fits", c.mask, header)

#makes a list of all images in the folder
images = glob.glob('pes*.fits')

#this loop makes a list of all the frames that are science frames (without sky frames)
sci_im = []
for i in images: 
    im  = pyfits.open(i)
    typ = im[0].header['IMAGETYP'] #This will be object for science frames but also for sky frames 
    obj = im[0].header['OBJECT'] #Use this to tell whether object frames are science or sky frames
    dith = string.split(obj,' ') #Split at the spaces. Science frames with have dith[0]=='dither'
    if dith[0] in standard_stars:
        print '\n'
        print 'skipping '+dith[0]+' '+i
        print '\n'
    else:
        if len(dith) > 1:
            if typ == 'object' and dith[1] == 'dither': #this gives science frames
                sci_im.append(i)

rmcosmicfits(sci_im)

