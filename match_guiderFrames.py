import numpy as np
import pyfits 
import string
import glob
import subprocess
import os.path as op

#This script assumes you have run sky subtraction on your science frames 

#********************#
# User defined paths #
#********************#

#path to science images (In reduction directory this is your data folder if old CURE, object folder if new CURE)
sci_path_base = '/Users/Briana/Documents/Grad_School/VIRUS_P/VP_MAR_2017/20170330/data'
#path to the guider files
guide_path = '/Users/Briana/Documents/Grad_School/VIRUS_P/VP_MAR_2017/20170330/guider'

#prefix for sky subtracted frames
frame_prefix = 'S*'

#+++++++++++++++++++ end user defined variables +++++++++++++++++++++#

#********************************#
# Find data and define functions #
#********************************#

#makes a list of all of the guider frames 
guide_im_list = glob.glob(op.join(guide_path,'*.fits'))

#makes a list of all of the science object names
object_list = glob.glob(op.join(sci_path_base,'*'))

#function that returns the date and time fromt header files 
def get_time(header):
    timeh = header['UT']
    dateh = header['DATE-OBS']

    t = string.split(timeh,':')
    hour = float(t[0])
    minu = float(t[1])
    secs = float(t[2]) 

    d = string.split(dateh,'-')
    year = int(d[0])
    month = int(d[1])
    day = int(d[2])

    time = hour + (minu/60) + (secs/3600) #in hours
    date = year+month+day

    return time, date

#*************************************#
# Build times lists for guider frames #
#*************************************#

gtime = []
gdate = []

#creates a list of dates and times for all of the guider frames
for i in range(len(guide_im_list)):
    g_im = pyfits.open(guide_im_list[i])
    hdr = g_im[0].header
 
    time = get_time(hdr)[0]
    date = get_time(hdr)[1]

    gtime.append(time)
    gdate.append(date)

np.save(op.join(guide_path,'guider_dates'),gdate)
np.save(op.join(guide_path,'guider_times'),gtime)
np.save(op.join(guide_path,'guider_names'),guide_im_list)

#*************************************#
# Match sci frames with guider frames #
#*************************************#

for o in range(len(object_list)):

    obj = op.basename(object_list[o])
    if (obj != 'trace') and (obj != 'flat'):
    
        print obj

        #path to science object
        sci_path = op.join(sci_path_base, obj)
        sci_im_list = glob.glob(op.join(sci_path, frame_prefix+'.fits'))

        file1 = 'sci_im_lis.lis'
        f1 = open(op.join(sci_path, file1),'w')
        for i in sci_im_list:
            f1.write(i + '\n') 
        f1.close()
        #subprocess.call('ls '+sci_path+'/Sp*.fits > '+sci_path+'/sci_im_lis.lis', shell=True)

        gtime = np.load(op.join(guide_path, 'guider_times.npy'))
        gdate = np.load(op.join(guide_path, 'guider_dates.npy'))
        gname = np.load(op.join(guide_path, 'guider_names.npy') )

        num_dithers = len(sci_im_list)

        stime = []
        sdate = [] 

        #get exposure time
        tim_im = pyfits.open(sci_im_list[0])
        tim_hdr = tim_im[0].header
        sci_exp = int(tim_hdr['EXPTIME']) #in seconds
        print 'Exposure Time: '+str(sci_exp/60)+' minutes'

        exp_length = sci_exp/float(3600) #in hours

        #creates a list of dates and times for all of the science frames
        for i in range(len(sci_im_list)):
            s_im = pyfits.open(sci_im_list[i])
            hdr = s_im[0].header

            hd_times = get_time(hdr)
         
            time = hd_times[0]
            date = hd_times[1]

            stime.append(time)
            sdate.append(date)

        for d in range(num_dithers):
            dmatch = np.where(gdate == sdate[d])[0]
            #tmatch = where(stime[d] < gtime[dmatch] < (stime[d]+exp_length))
            tmatch = np.where(np.logical_and(gtime[dmatch] > stime[d], gtime[dmatch] < (stime[d]+exp_length)))[0]

            im_names = gname[dmatch[tmatch]]

            print 'dither '+str(d)+': '+str(len(im_names))+' guider frames'

            filename = op.join(sci_path, 'time_'+str(d)+'.lis')
            print filename
            f2 = open(filename,'w')
            for i in range(len(im_names)):
                f2.write(op.join(guide_path, im_names[i])+'\n') 
            f2.close()

        print '\n'









