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
from virus_p_config import * 

#################################
# Defining which unit to reduce #
#################################

print ('#######################################')
print ('## RUNNING DATA REDUCTION ON VIRUS-P ##')
print ('#######################################')

###################
# Setting CUREBIN #
###################

CUREBIN = None
if not CUREBIN:
    CUREBIN = environ.get('CUREBIN')
if not CUREBIN:
    print("Please set CUREBIN as  environment variable or in the script")
    sys.exit(1)

###################################
# Defining which functions to run #
###################################

#if basic reduction is run need to specify specific routines to run 
# divide pixel flat and masterdark are not being used now
if basic:
    dividepf        = False
    normalize       = True 
    masterdark      = True
    masterarc       = True  
    mastertrace     = True 
else:
    dividepf        = False
    normalize       = False  
    masterdark      = False
    masterarc       = False  
    mastertrace     = False
    fixorange       = False

# This makes sure that the redux folder is only overwritten if the user chooses to run basic reduction
# If you user only wants to run deformer, skysubtract, fiberextract, or mkcube it used the data in redux 
if basic:
    all_copy = True
    RESTART_FROM_SCRATCH = True
else:
    all_copy = False
    RESTART_FROM_SCRATCH = False

#residual parameters from that do not need to be changed for automated use 
usemapping       = False # if headers are messed up, manual map IFUSLOT to SPECID 
initial_base     = ''

########################
# Specifying CURE opts #
########################

#specify opts for CURE routines used in this script
meanfitsopts    = "--new_error -k 3"
headfitsopts    = "-m -k EXPTIME -v 1 -V"
darkopts        = "--maxmem 1024 -s -t -m -k 2.8" 
arcopts         = "--maxmem 1024 -s -t -m -k 2.8"
traceopts       = "--maxmem 1024 -s -t -m -k 2.8"
deformeropts    = "-p 8 -C 8 -I [1.8,0.12,0.0,2.0] --debug --dump_psf_data"
subskyopts      = "-J --output-both -w "+str(window_size)+" -k "+str(sky_kappa[0])+','+str(sky_kappa[1])+" -m "+str(smoothing)+" -T "+str(sn_thresh)
fibextractopts  = "-P"
if diffAtmRef:
    cubeopts        = "-a "+str(sky_sampling)+" -k "+str(max_distance)+" -s "+str(cube_sigma)+" -d "
else:
    cubeopts        = "-a "+str(sky_sampling)+" -k "+str(max_distance)+" -s "+str(cube_sigma)

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
#DIR_DICT     = {    0:zro_dir,    1:drk_dir,    2:cmp_dir,    3:flt_dir,    4:sci_dir } # Dictionary for looping through directory names
DIR_DICT     = {    0:zro_dir,    1:cmp_dir,    2:flt_dir,    3:sci_dir }
#############################
# Define config directories #
#############################

#specifies directories for lines and mapping/cen files for LRS2 in the LRS2 config directory 
map_file = op.join(CUREBIN, '../config/IFUcen_VP2_27m.txt')

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
            self.clean                  = CLEAN_AFTER_DONE

            self.actionbase = initial_base  
            
            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            hdulist           = pyfits.open ( filename )
            self.trimsec      = "\"" + re.split('[\[ \] ]',hdulist[0].header['TRIMSEC'])[1] + "\""
            self.biassec      = "\"" + re.split('[\[ \] ]',hdulist[0].header['BIASSEC'])[1] + "\""

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


#Class for building dither file for CURE's mkcube
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
        s.append("#          _{L,R}.fits is added for the left and right spectrographs")
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
        return os.system(op.join(CUREBIN, command) +' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CUREBIN, command))
        
        
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
    

def extractfits(trimsec, frames):

    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]
    
    command = 'extractfits -r %s' % (trimsec)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command ( command, 0 )
    
    [f.addbase ('e') for f in frames] 
    
    return command


def addphotonnoise(frames):

    command = 'addphotonnoise'
    
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
    
    
def meanlampfits(dest_dir, basename , opts, frames):
    
    filenames = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in frames]

    mastername = op.join ( dest_dir , basename + '.fits')
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    run_cure_command( command, 0 )
    
    return command
 
def deformer_master(mastertrace,masterarc,linesfile,wave_range,ref_line,binning,opts):

    if binning == '2':
        frange = '0,1025'
    else:
        frange = '0,2065'

    wave_range = str(wave_range[0])+','+str(wave_range[1])
    
    command = 'deformer %s -s %s -W %s -F %s -o \"%s\" -l %s -a %s %s' % (opts,ref_line,wave_range,frange,redux_dir,linesfile,masterarc,mastertrace)  

    run_cure_command( command, 0 )

    return command   

def deformer(mastertrace,masterarc,linesfile,wave_range,ref_line,binning,opts):

    if binning == '2':
        frange = '0,1025'
    else:
        frange = '0,2065'

    wave_range = str(wave_range[0])+','+str(wave_range[1])
    
    command = 'deformer %s -s %s -W %s -F %s -o \"%s\" -l %s -a %s %s' % (opts,ref_line,wave_range,frange,redux_dir,linesfile,masterarc,mastertrace)  

    run_cure_command( command, 0 )

    return command
    
    
def subtractsky(frames,distmodel,fibermodel,opts,skymaster=""):
    
    filenames = [(redux_dir + '/' + sci_dir + '/pses' + f.basename + '.fits') for f in frames]

    for f in filenames:
        
        command = 'subtractsky %s %s -d %s -f %s %s' % (opts,skymaster,distmodel,fibermodel,f)  

        run_cure_command( command, 0 )

    return command

    
def fibextract_Resample(frames,base,distmodel,fibermodel,wave_range,nsample,use_ap_corr,opts):

    filenames = [(redux_dir + '/' + sci_dir + '/' + base + f.basename + '.fits') for f in frames]

    wave_range = str(wave_range[0])+','+str(wave_range[1])

    if use_ap_corr:
        ap_corr = '-c'
    else:
        ap_corr = ''

    for f in filenames:
    
        command = 'fiberextract %s %s -p FeR -W %s -n %s -d %s -f %s %s' %(opts,ap_corr,wave_range,nsample,distmodel,fibermodel,f)
        
        run_cure_command( command, 0 )

    return command


def fibextract(frames,base,distmodel,fibermodel,use_ap_corr,opts):

    filenames = [(redux_dir + '/' + sci_dir + '/' + base + f.basename + '.fits') for f in frames]

    wave_range = str(wave_range[0])+','+str(wave_range[1])

    if use_ap_corr:
        ap_corr = '-c'
    else:
        ap_corr = ''

    for f in filenames:
    
        command = 'fiberextract %s %s -x -d %s -f %s %s' %(opts,ap_corr,distmodel,fibermodel,f)
        
        run_cure_command( command, 0 )

    return command

def mkcube(ifufile,ditherfile,outname,opts):
    
    command = 'mkcube %s -o %s -i %s %s' %(opts,outname,ifufile,ditherfile)
        
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
    if redux_dir is None:
        print ( "Please provide a reduction directory \"redux_dir\"." )
        return None
    else:
        if not op.exists ( redux_dir ):
            os.mkdir ( redux_dir )
        else:
            if RESTART_FROM_SCRATCH:
                shutil.rmtree ( redux_dir )
                os.mkdir ( redux_dir )
        
    if DIR_DICT is None:        
        print ( "Please provide a directory dictionary \"DIR_DICT\"." )
        print ( "The directory dictionary order should match the file location directory." )
        return None

    #make a list of all of the raw files in the date directory 
    filenames = glob.glob ( date_folder + "/*.fits" ) 

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
        if not op.exists ( op.join ( redux_dir, DIR_DICT[i] ) ):
            os.mkdir ( op.join ( redux_dir, DIR_DICT[i] ) )
        #Loop through the files and find the ones that match the DIR_DICT data type
        typ_files = [f for f in filenames if file_typ_dict[f] == DIR_DICT[i]]  
        print (DIR_DICT[i]+': Found '+str(len(typ_files))+' files')

        # Loop through the retrieved files names to copy to new structure
        # Create a VirusFrame class for each frame that can maintain info for each original frame
        # The VirusFrame is critical for the rest of the reduction pipeline
        if all_copy:
            for t in typ_files:         
                shutil.copy ( t, op.join ( redux_dir, DIR_DICT[i] ) )
            for t in typ_files:        
                a = VirusFrame( op.join( redux_dir, DIR_DICT[i], op.basename ( t ) ) ) 
                vframes.append(copy.deepcopy(a))
        else:
            for t in typ_files:            
                a = VirusFrame(  t  ) 
                vframes.append(copy.deepcopy(a))
                        
    return vframes 

def basicred(redux_dir, DIR_DICT, basic = False, dividepf = False,
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
    9) divide pixelflat from cmp/flt/sci frames
    10) normalize cmps and flts 
    11) combine cmps and flts into masterarc and mastertrac
    '''

    print ('*************************')
    print ('* BUILDING IMAGE FRAMES *')
    print ('*************************')

    vframes = initial_setup (redux_dir, DIR_DICT )
    oframes = [v for v in vframes if v.type != "zero"] # gives "flt", "drk", "cmp", and "sci" frames (basically just not "zro")
    dframes = [v for v in vframes if v.type == "dark"] # gives dark frames
    cframes = [v for v in vframes if v.type == "flat" or v.type == "comp"] # gives "flt" and "cmp" frames
    lframes = [v for v in vframes if v.type == "comp"] # gives just "cmp" frames
    fframes = [v for v in vframes if v.type == "flat"] # gives just "flt" frames
    sframes = [v for v in vframes if v.type == "object"] # gives just "sci" frames

    file_binning = vframes[0].binning[0] #This will equal 2 if the files are binned, assumed to be the same for all files

    #make a copy of the virus_p_config file to be added to your directory
    #if the file exists - remove the file and replace it.
    if os.path.isfile(redux_dir+'/virus_p_config_'+redux_dir+'_copy.py') == True:
        os.remove(redux_dir+'/virus_p_config_'+redux_dir+'_copy.py')
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/virus_p_config.py', redux_dir+'/virus_p_config_'+redux_dir+'_copy.py' )
    else:
        shutil.copy ( os.path.dirname(os.path.realpath(__file__))+'/virus_p_config.py', redux_dir+'/virus_p_config_'+redux_dir+'_copy.py' )


    if basic:
        trimsec = vframes[0].trimsec # Trimsec assumed to be the same for all frames 
        biassec = vframes[0].biassec # Biassec assumed to be the same for all frames
        print ('************************')
        print ('* MAKE ERROR FRAME FOR *')
        print ('************************')
        mkerrorframe ( vframes)               # for all frames
        print ('*************************')
        print ('* SUBTRACT OVERSCAN FOR *')
        print ('*************************')
        subtractoverscan( biassec, vframes )   # for all frames
        print ('***************************')
        print ('* EXTRACT DATA REGION FOR *')
        print ('***************************')
        extractfits ( trimsec, vframes)       # for all frames
        
        vframesselect  = [v for v in vframes if v.type == "zero"] 
        print ('************************')
        print ('* BUILD MASTERBIAS FOR *')
        print ('************************')
        meanframefits   ( redux_dir, 'masterbias', meanfitsopts, vframesselect ) # meanfits for masterbias for unique specid
        masterbiasname = op.join ( redux_dir, 'masterbias.fits' ) 
        #oframesselect  = [o for o in oframes] 
        print ('***************************')
        print ('* SUBTRACT MASTERBIAS FOR *')
        print ('***************************')
        subtractbias   ( oframes, masterbiasname) # for all frames
    
        print ('**********************')
        print ('* MULTIPYING BY GAIN *')
        print ('**********************')
        addphotonnoise ( oframes ) # for all frames

    # Dividing by Pixel Flat
    if dividepf:
        print ('*************************')
        print ('* DIVDING BY PIXEL FLAT *')
        print ('*************************')
        pflat = op.join( pixflatdir, "pixelflat_cam{:03d}_{:s}.fits".format( uca, side.upper() ) )
        opt = "--file {:s}".format(pflat)
        dividepixelflat(oframes)
    
    # Create Master Dark Frames
    if masterdark:
        print ('***********************')
        print ('* BUILDING MASTERDARK *')
        print ('***********************')
        #dframesselect = [d for d in dframes if d.specid == uca] 
        meanframefits ( redux_dir, 'masterdark', meanfitsopts, dframes ) # meanfits for masterdark for unique specid

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
            meanlampfits(redux_dir, 'masterarc' , arcopts, lframes) 
        else:
            print ('Only one arc image: copying to masterarc')
            filename = [op.join ( f.origloc, f.actionbase + f.basename +'.fits') for f in lframes]
            mastername = op.join ( redux_dir , 'masterarc.fits')
            shutil.copy ( filename[0], mastername )
            efilename = [op.join ( f.origloc, 'e.' + f.actionbase + f.basename + '.fits') for f in lframes]
            emastername = op.join ( redux_dir , 'e.masterarc.fits')
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
            meanlampfits(redux_dir, 'mastertrace', traceopts, fframes)
        else:
            print ('Only one flat image: copying to mastertrace')
            filename = [op.join ( f.origloc, f.actionbase + f.basename + '.fits') for f in fframes]
            mastername = op.join ( redux_dir , 'mastertrace.fits')
            shutil.copy ( filename[0], mastername )
            efilename = [op.join ( f.origloc, 'e.' + f.actionbase + f.basename + '.fits') for f in fframes]
            emastername = op.join ( redux_dir , 'e.mastertrace.fits')
            shutil.copy ( efilename[0], emastername )

        for f in fframes:
            if f.clean:
                filename   = op.join ( f.origloc,        f.actionbase + f.basename + '.fits')
                filename_e = op.join ( f.origloc, 'e.' + f.actionbase + f.basename + '.fits')
                os.remove ( filename )
                os.remove ( filename_e )

    # Run Deformer
    if run_deformer:
        print ('*************************************************************************')
        print ('* RUNNING DEFORMER TO BUILD DISTORTION SOLUTION AND WAVELENGTH SOLUTION *')
        print ('*************************************************************************')
        shutil.copy ( lines_file, op.join(redux_dir,op.basename(lines_file)) )
        mastertrace = op.join ( redux_dir, 'mastertrace.fits' )
        masterarc   = op.join ( redux_dir, 'masterarc.fits' )
        linefile    = op.join ( redux_dir, op.basename(lines_file) )
        deformer ( mastertrace, masterarc, linefile, spec_range, ref_line, file_binning, deformeropts)
    
    # Run sky subtraction            
    if subsky:  
        print ('************************************************')
        print ('* PERFORMING SKY SUBTRACTION ON SCIENCE FRAMES *')
        print ('************************************************')
        distmodel = op.join ( redux_dir, 'mastertrace.dist' )
        fibermodel = op.join ( redux_dir, 'mastertrace.fmod' )

        print ('    ++++++++++++++++++++++++++++++')
        print ('     Sky Subtract Standard Frames ')
        print ('    ++++++++++++++++++++++++++++++')
        sframesselect = [s for s in sframes if s.object in standard_frames]
        print ('Found '+str(len(sframesselect))+' Standard Frames')
        subtractsky(sframesselect,distmodel,fibermodel,subskyopts)

        print ('    ++++++++++++++++++++++++++++++')
        print ('     Sky Subtract Science Frames ')
        print ('    ++++++++++++++++++++++++++++++') 
        #This is will equal 0 if there are no sky frames 
        useskyframe = len(sky_frames)
        #if there are no sky frames it just does normal sky subtraction 
        if useskyframe == 0:
            print ('No sky frames found: Fibers in frame will be used for sky model')
            sframesselect = [s for s in sframes if s.object in science_frames]
            print ('Found '+str(len(sframesselect))+' Science Frames')
            subtractsky(sframesselect,distmodel,fibermodel,subskyopts)
        #if there are sky frames 
        else:
            print ('Sky frames will be used for sky model')
            sky_ims = [s for s in sframes if s.object in sky_frames]
            sky_times = [t.datetime for t in sky_ims]
            sframesselect = [s for s in sframes if s.object in science_frames]
            print ('Found '+str(len(sframesselect))+' Science Frames')
            for f in sframesselect:
                sci_date = f.datetime
                closest_index = min(range(len(sky_times)), key=lambda i: abs(sky_times[i]-sci_date))
                skyframe = sky_ims[closest_index]
                skyframe_name = redux_dir + '/' + sci_dir + '/pses' + skyframe.basename + '.fits'

                skymaster = '-X '+skyframe_name+' -D '+distmodel+' -F '+fibermodel
                subtractsky(sframesselect,distmodel,fibermodel,subskyopts,skymaster)


    # Run fiberextract
    if fiberextract:  
        print ('****************************************')
        print ('* EXTRACTING SPECTRA IN SCIENCE FRAMES *')
        print ('****************************************')

        #finds if there are sky subtracted files. If so it uses those.
        Sfiles = glob.glob(redux_dir + "/" + sci_dir + "/Sp*.fits")
        if len(Sfiles) == 0:
            base = 'pses'
        else:
            base = 'Spses'

        science_objects = science_frames + standard_frames
        sframesselect = [s for s in sframes if s.object in science_objects]
        print ('Found '+str(len(sframesselect))+' Science Frames for Fiber Extraction')

        distmodel = redux_dir + "/mastertrace.dist"
        fibermodel = redux_dir + "/mastertrace.fmod"

        if wl_resample:
            print ('    ++++++++++++++++++++++++++')
            print ('     Resampling in Wavelength ')
            print ('    ++++++++++++++++++++++++++')

            if file_binning == '2':
                nsample = 1024
            else:
                nsample = 2048

            fibextract_Resample(sframesselect,base,distmodel,fibermodel,spec_range,nsample,use_ap_corr,fibextractopts) 
        else:
            print ('    +++++++++++++++++++++++++++++++')
            print ('     Extraction Without Resampling ')
            print ('    +++++++++++++++++++++++++++++++')

            fibextract(sframesselect,base,distmodel,fibermodel,use_ap_corr,fibextractopts)

    #CURE saves these files from deformer outside of the redux directory for some reason.
    #This moves them inside of the redux directory.
    left_files = glob.glob('*.log') + glob.glob('*.residuals')
    if len(left_files) > 0:
        for l in left_files:
            os.rename(l, op.join(redux_dir,l))

    #Run mkcube
    if makecube:
        print ('***********************')
        print ('* BUILDING DATA CUBES *')
        print ('***********************')
        location_prefix = redux_dir + "/" + sci_dir + "/" 
        os.chdir(location_prefix)

        #builds a list of files to build a data cube 
        #checks for if there are wavelength resampled fiber extracted and sky subtracted files
        Fefiles = glob.glob("FeRS*_sci_*.fits")
        #if not then checks for fiber extraced and sky subtracted files 
        if len(Fefiles) == 0:
            Fefiles = glob.glob("FeS*_sci_*.fits")
        #if not checks for any type of fiber extracted files (resampled or not)
        if len(Fefiles) == 0:
            Fefiles = glob.glob("Fe*_sci_*.fits")

        ditherfile = 'dither_vp.txt'

        for f in Fefiles:
            im  = pyfits.open(f)
            hdr = im[0].header
            #extracting header information for to know what channel file corresponds to and for dither file information 
            airmass = hdr['AIRMASS']


            #fix becuase deformer won't run on blue channel but need both sides for mkcube 
            if os.path.isfile(f[0:-7]+'_L.fits') == False:
                shutil.copy(f,f[0:-7]+'_L.fits')
                shutil.copy('e.'+f,'e.'+f[0:-7]+'_L.fits')


            psf      = 1.5
            basename = f[2:-7]
            outname  = f[0:-5]

            ditherf = open(ditherfile, 'w')
            ditherinfo.writeHeader(ditherf)
            ditherinfo.writeDither(ditherf,basename,"../mastertrace_"+str(uca),0.0,0.0,psf,1.00,airmass)

            mkcube(IFUfile,ditherfile,outname,cubeopts)  

            #fix becuase deformer won't run on blue channel but need both sides for mkcube
            if (uca == 501) and (side == 'R'):
                 os.remove(f[0:-7]+'_L.fits') 
                 os.remove('e.'+f[0:-7]+'_L.fits')  

    return vframes
    
def main():
    frames = basicred( redux_dir, DIR_DICT, basic = basic, dividepf = dividepf,
                      normalize = normalize, masterdark = masterdark, masterarc = masterarc, mastertrace = mastertrace )                 
    
if __name__ == '__main__':
    main()  
