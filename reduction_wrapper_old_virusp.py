# -*- coding: utf-8 -*-
"""
Created on Friday April  8 02:10:02 2016

This is a work in progress ... 

@author: gregz, brianaindahl
"""
from __future__ import print_function

import numpy as np
from scipy import interpolate
import pyfits
import glob
import six
import datetime as dt
import shutil
import sys
import os
import copy
import os.path as op
from os import environ
import re
import string
import cosmics
import argparse as ap
import importlib
#from virus_p_config import * 

#############################
# Define config file to use #
#############################

def parse_args(argv=None):
    """Parse the command line arguments
    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used
    Returns
    -------
    Namespace
        parsed arguments
    """
    description = "config_file"
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.ArgumentDefaultsHelpFormatter)
                        
    parser.add_argument("-config", nargs='?', type=str, help='''config file. "path_to_config/lrs2_config.py"''', default="virus_p_config.py")

    args = parser.parse_args(args=argv) 
           
    return args

args = parse_args()

config_file_name = args.config
config_arg       = config_file_name.split(".")[0]

config = importlib.import_module(config_arg, package=None)

#################################
# Defining which unit to reduce #
#################################

print ('#######################################')
print ('## RUNNING DATA REDUCTION ON VIRUS-P ##')
print ('#######################################')

###################
# Setting CUREVP #
###################

#Old CURE for VIRUS-P

CUREVP = None
if not CUREVP:
    CUREVP = environ.get('CUREVP')
if not CUREVP:
    print("Please set CUREVP as  environment variable or in the script")
    sys.exit(1)

#New CURE for VIRUS-P

CUREVP = None
if not CUREVP:
    CUREVP = environ.get('CUREVP')
if not CUREVP:
    print("Please set CUREVP as  environment variable or in the script")
    sys.exit(1)

#checking that LRS2 is defined in specconf.h 
cureversion = os.popen(op.join(CUREVP, 'cureversion')).readlines()
spec_define = cureversion[4].split(' ')[1]
instrument = spec_define.rstrip('\n')

if instrument == 'VIRUS-P':
    print ('CURE is set for VIRUS-P reduction')
else: 
    print ('You need to update specconf.h in CURE to define VIRUS-P')
    sys.exit('Right now specconf.h defines '+instrument)

###################################
# Defining which functions to run #
###################################

#if basic reduction is run need to specify specific routines to run 
# divide pixel flat and masterdark are not being used now
if config.basic:
    rmcosmics       = config.rmCosmics 
    normalize       = True 
    masterdark      = True
    masterarc       = True  
    mastertrace     = True 
    sort_sci        = True

else:
    rmcosmics       = False
    normalize       = False  
    masterdark      = False
    masterarc       = False  
    mastertrace     = False
    sort_sci        = False

# This makes sure that the redux folder is only overwritten if the user chooses to run basic reduction
# If you user only wants to run deformer, skysubtract, fiberextract, or mkcube it used the data in redux 
if config.basic:
    all_copy = True
    RESTART_FROM_SCRATCH = True
else:
    all_copy = False
    RESTART_FROM_SCRATCH = False

########################
# Specifying CURE opts #
########################

#specify opts for CURE routines used in this script
meanfitsopts    = "--new_error -k 3"
headfitsopts    = "-m -k EXPTIME -v 1 -V"
darkopts        = "--maxmem 1024 -s -t -m -k 2.8" 
arcopts         = "--maxmem 1024 -s -t -m -k 2.8"
traceopts       = "-r 400:550,1042:1144 -k 2.8 --norm_kappa 10"
deformeropts    = "--ref-files1 --ref-files2"
subskyopts      = "--output-both -w "+str(config.window_size)+" -k "+str(config.sky_kappa[0])+','+str(config.sky_kappa[1])+" -m "+str(config.smoothing)+" -T "+str(config.sn_thresh)
fibextractopts  = "-P"
cubeopts        = "-a "+str(config.airmass_lis)+" -k "+str(config.max_distance)+" -s "+str(config.cube_sigma)

#########################
# Defining data folders #
#########################

#specify folders where data types are stored
zro_dir  = "zero"
flt_dir  = "flat"
sci_dir  = "object"
cmp_dir  = "comp"
drk_dir  = "dark"

##########################
# Building Spec Libaries #
##########################

#dictionary of data type folders 
DIR_DICT     = {    0:zro_dir,    1:drk_dir,    2:cmp_dir,    3:flt_dir,    4:sci_dir } # Dictionary for looping through directory names

#############################
# Define config directories #
#############################

#specifies directories for lines and mapping/cen files for LRS2 in the LRS2 config directory 
map_file = op.join(CUREVP, '../config/IFUcen_VP2_27m.txt')

##################
# Define Classes #
##################

#class to build VIRUS frames for all data 
class VirusFrame:
    def __init__ ( self, filename = None):
        '''
        Initializing a VirusFrame for a given filename.
        This includes reading the header and setting reduction keywords.
        From this, the file with have attributes that can be tracked.
        '''
        
        if filename is None:
            print ( "No filename was given for VirusFrame to be initialized." )
            return None
        else:
            ######## OPEN AND KEEP FILES IN SOME SMART WAY ########
            self.filename               = filename
            self.origloc                = op.dirname  ( self.filename )
            self.basename, temp1        = op.basename ( self.filename ).split('.')
            self.clean                  = config.CLEAN_AFTER_DONE

            self.actionbase = ''  
            
            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            hdulist           = pyfits.open ( filename )
            #self.trimsec      = "\"" + re.split('[\[ \] ]',hdulist[0].header['TRIMSEC'])[1] + "\""
            #self.biassec      = "\"" + re.split('[\[ \] ]',hdulist[0].header['BIASSEC'])[1] + "\""

            obj_split        = hdulist[0].header['OBJECT'].split(' ')
            self.isdith      = len(obj_split)
            if self.isdith > 1:
                self.object      = obj_split[0]
                self.dither      = obj_split[2]
            else:
                self.object      = hdulist[0].header['OBJECT']
                self.dither      = 0

            self.type        = hdulist[0].header['IMAGETYP']            
            self.orggain     = hdulist[0].header['GAIN1']
            self.orgrdnoise  = hdulist[0].header['RDNOISE1']
            self.exptime     = hdulist[0].header['EXPTIME']
            self.binning     = hdulist[0].header['CCDSUM'].split(' ') 

            self.date        = hdulist[0].header['DATE-OBS']
            self.yr          = self.date.split('-')[0]
            self.mon         = self.date.split('-')[1]
            self.day         = self.date.split('-')[2]

            self.time        = hdulist[0].header['UT']
            self.hr          = self.time.split(':')[0]
            self.min         = self.time.split(':')[1]
            self.sec         = (self.time.split(':')[2]).split('.')[0]

            self.datetime    = dt.datetime(int(self.yr), int(self.mon), int(self.day), int(self.hr), int(self.min), int(self.sec))

    def addbase(self, action):
        if self.clean:
            filename   = op.join ( self.origloc,        self.actionbase + self.basename + '.fits')
            filename_e = op.join ( self.origloc, 'e.' + self.actionbase + self.basename + '.fits')            
            self.actionbase = action + self.actionbase
            if op.exists  ( filename ):
                os.remove ( filename )
            if op.exists  ( filename_e ):
                os.remove ( filename_e )                          
        else:
            self.actionbase = action + self.actionbase


class ditherinfo(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []
        s.append("# basename is the base file name of the data files.")
        s.append("# modelbase is the base file name of the dist, fmod, pmode files corresponding to the data files")
        s.append("# $Id:$")
        s.append("#")
        s.append("# basename modelbase ditherx dithery seeing norm airmass")
        s.append("#")
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeDither(cls, f, filename, basename, ditherx, dithery, seeing, norm, airmass):
        """Write something to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = ("%s %s %4.2f %4.2f %4.2f %4.2f %4.2f" %
             (filename, basename, ditherx, dithery, seeing, norm, airmass))
        f.write(s)
        f.write("\n")
        f.flush()

####################
# Define Functions #
####################

def run_cure_command(command, suppress_output=0, debug=1):
    '''
       Run a cure command
       
       Parameter suppress_output is used to quiet down. (default)
       Return the value of os.system(command)
    '''
    # store info in log?

    if debug:
        print('Running: \'%s\'' % command)
    if not suppress_output:
        return os.system(op.join(CUREVP, command) +' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CUREVP, command))


def updateheader(frames):
    
    filenames = [op.join ( f.origloc, f.basename + '.fits') for f in frames]
    #assume the gain and readnoise is the same for all images because they are all taken with the same detector
    gain1 = f.orggain
    read1 = f.orgrdnoise

    command1 = 'headfits -a -t float -k GAIN -v ' + str(gain1) + ' -V "same as GAIN1 added for CURE"'
    command2 = 'headfits -a -t float -k RDNOISE -v ' + str(read1) + ' -V "same as READNOISE1 added for CURE"'

    for i in xrange(len(filenames)):
        command1 = command1 + ' ' + filenames[i]
        command2 = command2 + ' ' + filenames[i]        

    run_cure_command( command1, 0 )
    run_cure_command( command2, 0 )
        
    return command1

        
def mkerrorframe(frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    command = 'mkerrorframe' 
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    run_cure_command( command, 0 )
        
    return command
    
        
def subtractoverscan(biassec, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    command = 'subtractfits -s -a -k 2.8 -t -o %s -z' % (biassec)

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('s') for f in frames] 

    return command
    
    
def subtractbias(frames, masterbiasname):
        
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    command = 'subtractfits -f %s' % ( masterbiasname ) 
        
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('s') for f in frames] 

    return command
    

def meanframefits(dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('master') for f in frames]

    return command

def meandarkfits(dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('master') for f in frames]

    return command

def subtractdark(frames, masterdarkname, drk_scale):
        
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]
        
    #scale the masterdark to the exposure time of the science frames
    command_1 = 'multiplyfits -p d%s_ -c %s' % ( drk_scale,drk_scale ) 
    command_1 = command_1 + ' ' + masterdarkname    
    run_cure_command( command_1, 0 )

    #subtract the scaled masterdark from the science frames 
    command_2 = 'subtractfits -p '' -f %s' % ( masterdarkname ) 
    command_2 = command_2 + ' ' + filenames
    run_cure_command( command_2, 0 )

    return command_2
    

def extractfits(trimsec, frames):

    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]
    
    command = 'extractfits -r %s' % (trimsec)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command ( command, 0 )
    
    [f.addbase ('e') for f in frames] 
    
    return command


def addphotonnoise(frames):

    command = 'addphotonnoise --gain_key GAIN1'
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    [f.addbase ('p') for f in frames] 
    
    return command
 
    
def dividepixelflat(frames, opts):

    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames if f.specid is uid]

    for i in filenames:
        command = 'dividefits %s %s' % (opts,filenames[i])
        run_cure_command( command, 0 )
    
    [f.addbase ('d') for f in frames] 
    
    return command 
    
    
def dividefits(filename, opts):

    command = 'dividefits %s %s' % (opts,filename)
    
    run_cure_command( command, 0 )
        
    return command
    
    
def addfits(filename, opts):

    command = 'addfits %s %s' % (opts,filename)
    
    run_cure_command( command, 0 )
        
    return command 
    

def addkeys(filename):

    command = 'addkeys %s' % (filename)
    
    command = command + ' << EOFK\n EXPTIME=1\n EOFK\n'
    
    run_cure_command( command, 0 )
    
    return command 
    
def headfits(filename,opts):

    command = 'headfits %s %s' % (opts,filename)
        
    run_cure_command( command, 0 )
    
    return command 


def meantracefits(dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '.fits')
    
    command = 'flatfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command
    
def meanlampfits(dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command


def rmcosmicfits(frames):

    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]
    
    for i in xrange(len(filenames)):

        print (filenames[i])

        array, header = cosmics.fromfits(filenames[i])
        # array is a 2D numpy array

        im_gain = header['GAIN1']
        im_RN = header['RDNOISE1']

        # Build the object :
        c = cosmics.cosmicsimage(array, gain=im_gain, readnoise=im_RN, satlevel = 65535.0, sigclip = 7.0, sigfrac = 0.3, objlim = 7.0)
        # There are other options, check the manual...

        # Run the full artillery :
        c.run(maxiter = 4)

        # Write the cleaned image into a new FITS file, conserving the original header :
        cosmics.tofits(filenames[i], c.cleanarray, header)

        # If you want the mask, here it is :
        #cosmics.tofits("s20160909T093737.7_066LL_sci_mask.fits", c.mask, header)

 
def deformer_master(mastertrace,masterarc,linesfile,wave_range,ref_line,opts):

    wave_range = str(wave_range[0])+','+str(wave_range[1])
    
    command = 'deformer3 -Q %s -o \"%s\" -l %s -a %s %s' % (opts,config.redux_dir,linesfile,masterarc,mastertrace)  

    run_cure_command( command, 0 )

    return command   


def deformer(sciframes,masterdist,masterfmod,outpath,wave_range,ref_line,opts):

    filenames = [op.join(f.origloc, f.basename + '.fits') for f in sciframes]

    wave_range = str(wave_range[0])+','+str(wave_range[1])
    
    command = 'deformer3 -Q %s -d %s -f %s -o \"%s\"' % (opts,masterdist,masterfmod,outpath)  

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]

    run_cure_command( command, 0 )

    return command
    
    
def subtractsky(frames,opts):

    for f in frames:

        filename = op.join(f.origloc, f.basename + '.fits')
        distname = op.join(f.origloc, f.basename + '.dist')
        fmodname = op.join(f.origloc, f.basename + '.fmod')
        
        command = 'subtractsky %s -d %s -f %s %s' % (opts,distname,fmodname,filename)  

        run_cure_command( command, 0 )

    return command
    
    
def subtractsky_frame(frame,skyframe,skyscale,opts):

    filename = op.join(frame.origloc, frame.basename + '.fits')
    distname = op.join(frame.origloc, frame.basename + '.dist')
    fmodname = op.join(frame.origloc, frame.basename + '.fmod')

    filesky = op.join(skyframe.origloc, skyframe.basename + '.fits')
    distsky = op.join(skyframe.origloc, skyframe.basename + '.dist')
    fmodsky = op.join(skyframe.origloc, skyframe.basename + '.fmod')

    skyframeopts = '-X '+filesky+' -D '+distsky+' -F '+fmodsky
    
    command = 'subtractsky %s --x-sky-scaling %s %s -d %s -f %s %s' % (opts,skyscale,skyframeopts,distname,fmodname,filename)  

    run_cure_command( command, 0 )

    return command
    
def fibextract_Resample(filenames,skysub_files,wave_range,nsample,opts):

    wave_range = str(wave_range[0])+','+str(wave_range[1])

    for f in filenames:

        if skysub_files:
            dist = op.dirname( f ) + '/' + op.basename( f ).split('.')[0][1::] + '.dist'
            fmod = op.dirname( f ) + '/' + op.basename( f ).split('.')[0][1::]+ '.fmod'
        else:
            dist = op.dirname( f ) + '/' + op.basename( f ).split('.')[0] + '.dist'
            fmod = op.dirname( f ) + '/' + op.basename( f ).split('.')[0] + '.fmod'

        command = 'fiberextract %s -p FeR -W %s -n %s -d %s -f %s %s' %(opts,wave_range,nsample,dist,fmod,f)
        
        run_cure_command( command, 0 )

    return command


def fibextract(filenames,skysub_files,opts):

    wave_range = str(wave_range[0])+','+str(wave_range[1])

    for f in filenames:

        dist = op.dirname( f ) + '/' + op.basename( f ).split('.')[0] + '.dist'
        fmod = op.dirname( f ) + '/' + op.basename( f ).split('.')[0] + '.fmod'
    
        command = 'fiberextract %s -x -d %s -f %s %s' %(opts,dist,fmod,f)
        
        run_cure_command( command, 0 )

    return command

def mkcube(ifufile,ditherfile,outname,diffatmref,opts):

    #the default for the DAR correction on CURE is True so adding -d turns off DAR correction
    if diffatmref:
        dar = ""
    else:
        dar = "-d"
    
    command = 'mkcube %s %s -o %s -i %s %s' %(opts,dar,outname,ifufile,ditherfile)
        
    run_cure_command( command, 0 )

    return command
               

def initial_setup (redux_dir = None, DIR_DICT = None):
    '''
    Running the initial setup which includes:
    1) Building a standard reduction folder structure
    2) Copying files from specified location into the structure
    3) Creating class variable, VirusFrame, for each file which records
    keywords for other functions as well as time of the observation and basename
    '''

    vframes = [] # will fill this list with VirusFrame class objects for each image
    if config.redux_dir is None:
        print ( "Please provide a reduction directory \"redux_dir\"." )
        return None
    else:
        if not op.exists ( config.redux_dir ):
            os.mkdir ( config.redux_dir )
        else:
            if RESTART_FROM_SCRATCH:
                shutil.rmtree ( config.redux_dir )
                os.mkdir ( config.redux_dir )
        
    if DIR_DICT is None:        
        print ( "Please provide a directory dictionary \"DIR_DICT\"." )
        print ( "The directory dictionary order should match the file location directory." )
        return None

    #make a list of all of the raw files in the date directory 
    filenames = glob.glob ( config.date_folder + "/*.fits" ) 

    #iterate through the files and save the file type to a list
    file_typ_lis = []
    for f in filenames:
        hdulist  = pyfits.open ( f )
        ftype    = hdulist[0].header['IMAGETYP'] 
        file_typ_lis.append(ftype)

    #build a dictionary of the filenames and thier types for sorting
    file_typ_dict = dict(zip(filenames, file_typ_lis))
    
    # Loop through the file location directories     
    for i in xrange ( len ( DIR_DICT ) ):
        # If the reduction location exists, don't re-make the directory (also, any files in that directory remain)
        if not op.exists ( op.join ( config.redux_dir, DIR_DICT[i] ) ):
            os.mkdir ( op.join ( config.redux_dir, DIR_DICT[i] ) )
        #Loop through the files and find the ones that match the DIR_DICT data type
        typ_files = [f for f in filenames if file_typ_dict[f] == DIR_DICT[i]]  

        # Loop through the retrieved files names to copy to new structure
        # Create a VirusFrame class for each frame that can maintain info for each original frame
        # The VirusFrame is critical for the rest of the reduction pipeline
        if all_copy:
            for t in typ_files:         
                shutil.copy ( t, op.join ( config.redux_dir, DIR_DICT[i] ) )
            for t in typ_files:        
                a = VirusFrame( op.join( config.redux_dir, DIR_DICT[i], op.basename ( t ) ) ) 
                vframes.append(copy.deepcopy(a))
        else:
            for t in typ_files:            
                a = VirusFrame(  t  ) 
                vframes.append(copy.deepcopy(a))
                        
    return vframes 

def basicred(redux_dir, DIR_DICT, basic = False,
              normalize = False, masterdark = False, masterarc = False, mastertrace = False):
    '''
    Running the basic reduction which includes:
    1) Overscan subtract zero frames
    2) Trim zero frames
    3) Create master bias frame
    4) Create blank error frame with readnoise in it
    5) Overscan subtract and trim cmp/flt/sci frames
    6) Subtract master bias from cmp/flt/sci frames
    7) ccdcombine frames which puts units in e-
    8) add photon noise to combined frames
    10) normalize cmps and flts 
    11) combine cmps and flts into masterarc and mastertrac
    '''

    print ('*************************')
    print ('* BUILDING IMAGE FRAMES *')
    print ('*************************')

    vframes = initial_setup (config.redux_dir, DIR_DICT )
    oframes = [v for v in vframes if v.type != "zero"] # gives "flt", "drk", "cmp", and "sci" frames (basically just not "zro")
    dframes = [v for v in vframes if v.type == "dark"] # gives dark frames
    cframes = [v for v in vframes if v.type == "flat" or v.type == "comp"] # gives "flt" and "cmp" frames
    lframes = [v for v in vframes if v.type == "comp"] # gives just "cmp" frames
    fframes = [v for v in vframes if v.type == "flat"] # gives just "flt" frames
    zframes = [v for v in vframes if v.type == "zero"]
    sframes = [v for v in vframes if v.type == "object"] # gives just "sci" frames

    #check that there are files for all of the needed cal data 
    #Script ends if missing data 
    if len(zframes) == 0:
        sys.exit('No bias frames provided')
    else:
        print ('Found '+str(len(zframes))+' bias frames')

    if len(fframes) == 0:
        sys.exit('No flat frames provided')
    else:
        print ('Found '+str(len(fframes))+' flat frames')

    if len(lframes) == 0:
        sys.exit('No comp frames provided')   
    else:    
        print ('Found '+str(len(lframes))+' comp frames')

    #darks are optional so if there are non the script can continue
    #If darks found and user choose subDarks: a masterdark is created and subtracted from science frames
    print ('Found '+str(len(dframes))+' dark frames')
    if len(dframes) != 0:
        masterdark      = config.subDarks
        subtractdark    = config.subDarks
    else:
        masterdark      = False
        subtractdark    = False    

    #Checks that there are science frames for each standar,sky or science object 
    #If none found for one object: script ends and prints all sci objects found
    print ('Found '+str(len(sframes))+' science frames')
    sci_obj_names = [s.object for s in sframes]     #finds all science objects in data provided 
    sci_obj_names = list(set(sci_obj_names))    #compresses list so one of each object 
    object_list = config.standard_frames + config.sky_frames + config.science_frames #combine all objects to iterate through
    for o in object_list:
        stframes = [s for s in sframes if s.object == o]
        if len(stframes) == 0:
            print ('No '+o+' frames found')
            sys.exit('Science objects found: '+str(sci_obj_names))
        else:
            print ('    Found '+str(len(stframes))+' '+o+' frames')

    file_binning = vframes[0].binning[0] #This will equal 2 if the files are binned, assumed to be the same for all files

    #set trimsec and biassec based on if the data is binned or unbinned
    #can't trust header values 

    if file_binning == '2':
        print ('FILES ARE BINNED')
        biassec="1029:1056,2:2047"  # Biassec assumed to be the same for all frames
        trimsec="1:1024,1:2048"     # Trimsec assumed to be the same for all frames 
    else:
        print ('FILES ARE UNBINNED')
        biassec="2053:2080,2:2047"  # Biassec assumed to be the same for all frames
        trimsec="1:2048,1:2048"     # Trimsec assumed to be the same for all frames 

    #make a copy of the virus_p_config file to be added to your directory
    #if the file exists - remove the file and replace it.
    if os.path.isfile(config.redux_dir+'/virus_p_config_'+config.redux_dir+'_copy.py') == True:
        os.remove(config.redux_dir+'/virus_p_config_'+config.redux_dir+'_copy.py')
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/virus_p_config.py', config.redux_dir+'/virus_p_config_'+config.redux_dir+'_copy.py' )
    else:
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/virus_p_config.py', config.redux_dir+'/virus_p_config_'+config.redux_dir+'_copy.py' )


    if basic:
        print ('*******************************')
        print ('* UPDATE HEADER KEYS FOR CURE *')
        print ('*******************************')
        updateheader ( vframes)               # for all frames
        print ('********************')
        print ('* MAKE ERROR FRAME *')
        print ('********************')
        mkerrorframe ( vframes)               # for all frames
        print ('*********************')
        print ('* SUBTRACT OVERSCAN *')
        print ('*********************')
        subtractoverscan( biassec, vframes )   # for all frames
        print ('***********************')
        print ('* EXTRACT DATA REGION *')
        print ('***********************')
        extractfits ( trimsec, vframes)       # for all frames

        # Remove cosmic rays using L.A.cosmic
        if rmcosmics:
            print ('***********************************')
            print ('* REMOVE COSMIC RAYS (SCI IMAGES) *')
            print ('***********************************')
            rmcosmicfits ( sframes )       # for sci frames - because this is slow                    
    
        print ('********************')
        print ('* BUILD MASTERBIAS *')
        print ('********************')
        meanframefits   ( config.redux_dir, 'masterbias', meanfitsopts, zframes ) # meanfits for masterbias
        masterbiasname = op.join ( config.redux_dir, 'masterbias.fits' ) 
        #oframesselect  = [o for o in oframes] 
        print ('***********************')
        print ('* SUBTRACT MASTERBIAS *')
        print ('***********************')
        subtractbias   ( oframes, masterbiasname) # for all frames
    
        print ('**********************')
        print ('* MULTIPYING BY GAIN *')
        print ('**********************')
        addphotonnoise ( oframes ) # for all frames
    
    # Create Master Dark Frames
    if masterdark:
        print ('*******************************************************')
        print ('* BUILDING MASTERDARK and SUBTRACTING FROM SCI FRAMES *')
        print ('*******************************************************')
        sci_exptime = list(set([float(s.exptime) for s in sframes]))
        drk_exptime = list(set([float(d.exptime) for d in dframes]))
        for e in sci_exptime:
            close_drk = min(drk_exptime, key=lambda x:abs(x-e))
            dframesselect = [d for d in dframes if d.exptime == close_drk] 
            meandarkfits (config.redux_dir, 'masterdark_'+close_drk+'sec', meanfitsopts, dframesselect ) # meanfits for masterdark for unique specid
       
            #Subtracts Master Dark from Science Frames  
            scale_drk = e/close_drk
            sframesselect = [s for s in sframes if s.exptime == e] 
            masterdarkname = op.join ( config.redux_dir, 'masterdark_'+close_drk+'sec' + '.fits' )
            subtractdark ( sframesselect, masterdarkname, scale_drk) # for sci frames 

    if normalize:
        print ('***************************************************')
        print ('* NORMALIZING CALIBRATION IMAGES BY EXPOSURE TIME *')
        print ('***************************************************')
        # Normalizing Calibration Frames by Exposure Time
        for cal in cframes:
            opt     = '-c {:0.1f}'.format(cal.exptime)
            filename = op.join ( cal.origloc, cal.actionbase + cal.basename + '.fits' )
            dividefits ( filename, opt )
            cal.addbase ('d')
            filename = op.join ( cal.origloc, cal.actionbase + cal.basename + '.fits' )
            headfits ( filename , headfitsopts )
            
    if masterarc:        
        # Making Master Arc Frames
        print ('**********************')
        print ('* BUILDING MASTERARC *')
        print ('**********************')

        if len(lframes)>1:
            meanlampfits(config.redux_dir, 'masterarc' , arcopts, lframes) 
        else:
            print ('Only one arc image: copying to masterarc')
            filename = [op.join ( f.origloc, f.actionbase + f.basename +'.fits') for f in lframes]
            mastername = op.join ( config.redux_dir , 'masterarc.fits')
            shutil.copy ( filename[0], mastername )
            efilename = [op.join ( f.origloc, 'e.' + f.actionbase + f.basename + '.fits') for f in lframes]
            emastername = op.join ( config.redux_dir , 'e.masterarc.fits')
            shutil.copy ( efilename[0], emastername )

        for l in lframes:
            if l.clean:
                filename   = op.join ( l.origloc,        l.actionbase + l.basename + '.fits')
                filename_e = op.join ( l.origloc, 'e.' + l.actionbase + l.basename + '.fits')
                os.remove ( filename )
                os.remove ( filename_e )
         
    # Making Master Trace Frames
    if mastertrace:
        print ('************************')
        print ('* BUILDING MASTERTRACE *')
        print ('************************')
        if len ( fframes ) > 1:
            meantracefits(config.redux_dir, 'mastertrace', traceopts, fframes)
        else:
            print ('Only one flat image: copying to mastertrace')
            filename = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in fframes]
            mastername = op.join ( config.redux_dir , 'mastertrace.fits')
            shutil.copy ( filename[0], mastername )
            efilename = [op.join ( f.origloc, 'e.' + f.actionbase + f.basename + '.fits') for f in fframes]
            emastername = op.join ( config.redux_dir , 'e.mastertrace.fits')
            shutil.copy ( efilename[0], emastername )

        for f in fframes:
            if f.clean:
                filename   = op.join ( f.origloc,        f.actionbase + f.basename + '.fits')
                filename_e = op.join ( f.origloc, 'e.' + f.actionbase + f.basename + '.fits')
                os.remove ( filename )
                os.remove ( filename_e )           

    # Run Master Deformer
    if config.run_deformer:
        print ('********************************************************************************')
        print ('* RUNNING DEFORMER TO BUILD MASTER DISTORTION SOLUTION AND WAVELENGTH SOLUTION *')
        print ('********************************************************************************')
        #check that basic has been run 
        trace_files = glob.glob(op.join(config.redux_dir,'mastertrace*'))
        arc_files   = glob.glob(op.join(config.redux_dir,'masterarc*'))
        if len(trace_files) == 0 or len(arc_files) == 0:
            sys.exit("You must run basic reduction before you can run deformer")

        #Run deformer master to create the master distortion solution 
        shutil.copy ( config.lines_file, op.join(config.redux_dir,op.basename(config.lines_file)) )
        mastertrace = op.join ( config.redux_dir, 'mastertrace.fits' )
        masterarc   = op.join ( config.redux_dir, 'masterarc.fits' )
        linefile    = op.join ( config.redux_dir, op.basename(config.lines_file) )
        deformer_master ( mastertrace, masterarc, linefile, config.spec_range, config.ref_line, deformeropts)

    if sort_sci:
        print ('**************************************************')
        print ('* SORTING SCIENCE FRAMES INTO OBJECT DIRECTORIES *')
        print ('**************************************************')  
        #iterate through science,standard, and sky frames and build folder for each   
        for s in object_list:
            os.mkdir ( op.join( config.redux_dir, sci_dir, s ))
        #For each science frame copy it into correct object directory within science folder (along with error frame)
        for s in sframes:
            sname = op.join ( s.actionbase + s.basename + '.fits')
            shutil.move (op.join(s.origloc, sname ), op.join( s.origloc, s.object ) )
            shutil.move (op.join(s.origloc, 'e.'+sname ), op.join( s.origloc, s.object ) ) 

    print ('**********************************')
    print ('* BUILDING SCIENCE FRAME OBJECTS *')
    print ('**********************************') 
    #This builds a new VIRUS frame for all science objects with updated path information
    sorted_sframes = glob.glob( op.join( config.redux_dir, sci_dir, '*' , '*.fits' ) )
    sciframes = []
    for f in sorted_sframes:
        #make sure not to include error frames in the list of image frames 
        if op.basename(f)[0:2] != 'e.':
            a = VirusFrame( f ) 
            sciframes.append( copy.deepcopy( a ) )

    orig_sci = [s for s in sciframes if s.basename[0:4] == 'pses'] #reduced science frames 

    #Run deformer on science frames using master deformer solution 
    if config.run_deformer:
        print ('**************************************')
        print ('* RUNNING DEFORMER ON SCIENCE FRAMES *')
        print ('**************************************')
        #Run deformer to create distortion solution for each of the science frames 
        dist_master = op.join ( config.redux_dir, 'mastertrace.dist' )
        fmod_master = op.join ( config.redux_dir, 'mastertrace.fmod' )
        #iterate through all science objects and run deformer - must iterated so output files get saved in sorted sci directory structure
        for o in object_list:
            objectframeselect = [s for s in orig_sci if s.object == o]
            outpath = op.join( config.redux_dir, sci_dir, o )
            deformer ( objectframeselect, dist_master, fmod_master, outpath, config.spec_range, config.ref_line, deformeropts )

    # Run sky subtraction            
    if config.subsky:  
        print ('************************************************')
        print ('* PERFORMING SKY SUBTRACTION ON SCIENCE FRAMES *')
        print ('************************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(config.redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run sky subtraction")

        print ('    ++++++++++++++++++++++++++++++')
        print ('     Sky Subtract Standard Frames ')
        print ('    ++++++++++++++++++++++++++++++')
        sframesselect = [s for s in orig_sci if s.object in config.standard_frames]
        print ('Found '+str(len(sframesselect))+' Standard Frames')
        subtractsky(sframesselect,subskyopts)

        print ('    ++++++++++++++++++++++++++++++')
        print ('     Sky Subtract Science Frames ')
        print ('    ++++++++++++++++++++++++++++++') 
        #This is will equal 0 if there are no sky frames 
        useskyframe = len(config.sky_frames)
        #if there are no sky frames it just does normal sky subtraction 
        if useskyframe == 0:
            print ('No sky frames found: Fibers in frame will be used for sky model')
            sframesselect = [s for s in orig_sci if s.object in config.science_frames]
            print ('Found '+str(len(sframesselect))+' Science Frames')
            subtractsky(sframesselect,subskyopts)
        #if there are sky frames 
        else:
            print ('Sky frames will be used for sky model')
            #list of sky frames 
            sky_ims = [s for s in orig_sci if s.object in config.sky_frames]
            #list of datetimes for each sky frame for comparison
            sky_times = [t.datetime for t in sky_ims]
            #list of science frames
            sframesselect = [s for s in orig_sci if s.object in config.science_frames]
            print ('Found '+str(len(sframesselect))+' Science Frames')
            #For each science frame: find the sky frame taken closest to the time of your observation
            for f in sframesselect:
                sci_date = f.datetime
                closest_index = min(range(len(sky_times)), key=lambda i: abs(sky_times[i]-sci_date))
                skyframe = sky_ims[closest_index]
                #define the scale factor to scale up sky exposure to to deal with different expsosure times between sky and sci frames
                #scale factor for the sky frame is the sci exposure time divided by they sky frame exposure time
                #this scales by exposure time and additionally by a factor the user chooses 
                skyscale = float(f.exptime)/float(skyframe.exptime) * config.sky_scaling
                subtractsky_frame(f,skyframe,skyscale,subskyopts)

    # Run fiberextract
    if config.fiberextract:  
        print ('****************************************')
        print ('* EXTRACTING SPECTRA IN SCIENCE FRAMES *')
        print ('****************************************')
        #check that deformer has been run 
        dist_files = glob.glob(op.join(config.redux_dir,'*.dist'))
        if len(dist_files) == 0:
            sys.exit("You must run deformer before you can run fiber extract")

        subsky_sci = glob.glob(config.redux_dir + "/" + sci_dir + "/*/" + "Spses*.fits")
        if len(subsky_sci) == 0:
            skysub_files = False
            sci_objs = config.science_frames + config.standard_frames
            sci_frames = [s for s in orig_sci if s.object in sci_objs]
            subsky_sci = [(f.origloc + '/' + f.basename + '.fits') for f in sci_frames]
        else:
            skysub_files = True

        print ('Found '+str(len(subsky_sci))+' Science Frames for Fiber Extraction')

        if config.wl_resample:
            print ('    ++++++++++++++++++++++++++')
            print ('     Resampling in Wavelength ')
            print ('    ++++++++++++++++++++++++++')

            if file_binning == '2':
                nsample = 1024
            else:
                nsample = 2048

            fibextract_Resample(subsky_sci,skysub_files, config.spec_range,nsample,fibextractopts) 
        else:
            print ('    +++++++++++++++++++++++++++++++')
            print ('     Extraction Without Resampling ')
            print ('    +++++++++++++++++++++++++++++++')

            fibextract(subsky_sci, skysub_files, fibextractopts)

    #CURE saves these files from deformer outside of the redux directory for some reason.
    #This moves them inside of the redux directory.
    left_files = glob.glob('*.log') + glob.glob('*.residuals')
    if len(left_files) > 0:
        for l in left_files:
            os.rename(l, op.join(config.redux_dir,l))

    #mkcube and collapse cube have NOT been edited to work 

    #Run mkcube
    if config.makecube:
        print ('***********************')
        print ('* BUILDING DATA CUBES *')
        print ('***********************')

        for s in config.science_frames:
        #cd inside of the science directory 
            location_prefix = op.join(config.redux_dir, sci_dir, s)
            os.chdir(location_prefix)
            Fefiles = glob.glob(op.join(config.redux_dir, sci_dir, s, "Fe*.fits"))

            if len(Fefiles) == 0:
                sys.exit('You must run fiber extraction before you can build data cubes')

            ditherfile = 'dither_vp.txt'

            for f in Fefiles:
                #    ditherx = [0.0,0.0,0.0,-1.197]
                #    dithery = [0.0,0.0,0.0,0.565]
                #    psf     = [2.0,2.0,2.0,2.0]
                #    norm    = [1.0,1.0,1.0,1.0]
                #    airmass = [1.23,1.23,1.23,1.23]
                airmass  = f.airmass
                psf      = 4.0
                basename = f.basename[2:-1] #????
                outname  = f.basename

                ditherf = open(ditherfile, 'w')
                ditherinfo.writeHeader(ditherf)
                ditherinfo.writeDither(ditherf,basename,"../mastertrace_"+str(uca),0.0,0.0,psf,1.00,airmass)

                mkcube(IFUfile,ditherfile,outname,config.diffAtmRef,cubeopts)    

            #cd back into the reduction directory 
            os.chdir('../../../')

    # Run collapse cube
    if config.collapseCube:
        print ('***************************************')
        print ('* COLLAPSING DATA CUBE TO BUILD IMAGE *')
        print ('***************************************')

        cube_sci = glob.glob(config.redux_dir + "/" + sci_dir + "/*/" + "Cu*.fits")

        for s in sci_objects:
            #cd into the science directory 
            location_prefix = config.redux_dir + "/" + sci_dir + "/" + s + "/" 

            #makes sure there are actually data cubes made from wavelength resampled, fiber extracted data in the sci directory 
            #If data cubes were made from fiber extracted fibers that do not wl resample they do not contain WCS info needed
            if len(Cufiles) == 0:
                sys.exit("You must build data cubes from wavelength resampled, fiber extracted data before running collapse cube")

            #user defined wavelength range to collapse cube 
            low_wave  = config.col_wave_range[0]
            high_wave = config.col_wave_range[1]

            #track number of cubes used in order to inform user if their values fall out of bounds and no images made
            num_cubes = 0

            min_wave_set = []
            max_wave_set = []
            #iterate through cube files 
            for c in Cufiles:
                im  = pyfits.open(c)
                hdr = im[0].header
                dat = im[0].data

                #read header for wavelength solution information
                lenZ  = np.shape(dat)[0]
                CRVAL = hdr['CRVAL3']                                        
                CDELT = hdr['CDELT3']
                Side  = hdr['CCDPOS']

                #build wavelength solution and find min and max wavelength of that solution
                wave_sol = np.add(np.arange(0,(lenZ*CDELT)+1,CDELT),CRVAL)
                max_wave = np.amax(wave_sol)
                min_wave = np.amin(wave_sol)

                #append the values for each frame to inform user of bounds of this data set if their values are out of bounds
                max_wave_set.append(max_wave)
                min_wave_set.append(min_wave)

                #If they choose to collapse entire cube ([0,0]) all data cubes are collapsed into images
                if low_wave == 0 and high_wave == 0:
                    num_cubes = num_cubes + 1 
                    print('Building image from '+spec_chan+' channel cube: '+op.basename(c))
                    print('Collapsing entire cube')

                    sum_image  = np.sum(dat, axis=0) #sums data cube in z direction
                    pyfits.writeto( location_prefix+'Col'+op.basename(c), sum_image, header=hdr, clobber=True)
                    print('\n')

                #If the user choose a wavelength range check if it in range of this cube 
                if (low_wave >= min_wave) and (high_wave <= max_wave):
                    num_cubes = num_cubes + 1 
                    print('Building image from cube: '+op.basename(c))
                    print('Collapsing cube from '+str(low_wave)+' to '+str(high_wave))

                    #find what index these wavelengths most closely correspond to. 
                    low_ind  = (np.abs(wave_sol-low_wave)).argmin()
                    high_ind = (np.abs(wave_sol-high_wave)).argmin()

                    #find wavelength value at this index - not exactly users choice so want to print value
                    low_val  = str(wave_sol[low_ind]).split('.')[0]
                    high_val = str(wave_sol[high_ind]).split('.')[0]

                    #Build image from that data cube 
                    dat_region = dat[low_ind:high_ind,:,:] #build region from low to high z - include all x,y
                    sum_image  = np.sum(dat_region, axis=0) #sum the image in the z direction
                    pyfits.writeto( location_prefix+'Col_'+low_val+'_'+high_val+'_'+op.basename(c), sum_image, header=hdr, clobber=True)
                    print('\n')

            #if num_cubes is 0: all cubes out of wavelength range of users choice 
            if num_cubes == 0:
                print ("Wavelength range you choose for collapse cube is out of range")
                sys.exit("This VIRUS-P data set ranges between "+str(np.amin(min_wave_set))+" and "+str(np.amax(max_wave_set))+" Angstroms")


    return vframes
    
def main():
    frames = basicred( config.redux_dir, DIR_DICT, basic = config.basic,
                      normalize = normalize, masterdark = masterdark, masterarc = masterarc, mastertrace = mastertrace )                 
    
if __name__ == '__main__':
    main()  
