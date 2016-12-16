############################################
# Configuration file for VIRUS-P reduction # 
############################################

#reduction_wrapper_virusP.py should be run inside of a folder containing the date folder containing raw data

redux_dir        = "VP_apr15_B1"	 	#name of the folder reduction is run 
date_folder      = "20150417"		#date folder containing raw data 

#Choose the spectral range of your VIRUS-P data 
spec_range   	 = [3500,5850] 		#Spectral range of that nights observations in angstroms

lines_file       = "/Users/Briana/Documents/Grad_School/VIRUS_P/apr15/lines.apr15"	#path to the lines file for this data set
IFU_file 		 = "/Users/Briana/Documents/cure_v2/cure/config/IFUcen_VP2_27m.txt" #path to the IFUcen file for this data set

CLEAN_AFTER_DONE = True 

########################
# Header Objects Names #
########################
#These should be lists of the objects names as written in the headers for the keyword OBJECT 
#Do NOT include 'dither #' after the object name even if it is part of the keyword OBJECT 
#standard and sky frames can be left as empty arrays if you don't have those files

#This should be a list of the objects names of standard star objects as written in the headers
standard_frames  = ['Feige34',]
#If you took sky frames, this should be a list of the object names of those frames as written in the headers
sky_frames		 = ['M82_sky',]
#This should be the object names of the science frames as written in the headers
science_frames      = ['M82_F1',]	

#####################
# REDUCTION OPTIONS #
#####################
basic           = False		#run basic reduction (normalize, build mastertrace and masterarc)
run_deformer    = False		#run deformer to map spectral traces and build wavelength solution for fiber extraction
subsky          = True  	#run sky subtraction on sci images - Need to have run deformer, only runs on non-extracted spectra
fiberextract    = False	 	#extract spectra and save into fits file - Need to have run deformer
makecube 		= False		#builds data cube out of fiber extracted image - Need to have run fiberextract
collapseCube 	= False		#collapse data cube to make an image of a wavelength range of the users choice

#------------#
# basic opts #
#------------#
ref_line		= 5 			#[integer] arc line in line list that is the brightest but not blended (used as reference line)
subDarks 		= False			#[True/False] If True darks will be subtracted from science images (default: False)
rmCosmics	 	= True  		#[True/False] If True the program L.A.Cosmic is used to eliminate cosmic rays (default: True)

#----------------------#
# sky subtraction opts #
#----------------------#
window_size 	= 250 			#[integer] Size (in image lines) of moving window for sky median
sky_kappa 		=[3.5,3.5]		#[floatarray] Lower and upper kappa for final sky kappa-sigma clipping.
smoothing 		= 5.0 			#[float] Smoothing factor for approximating spline.
sn_thresh		= 5 			#[float] Minimum signal to noise to flag fiber as continuum and ignore it during sky generation.

#--------------------#
# fiber extract opts #
#--------------------#
wl_resample 	 = True 		#[True/False] If True it will resample in wavelength, Does not resample in wavelength if False 
use_ap_corr	 	 = False		# Apply aperture correction using the fiber model (slow).

#----------------#
# make cube opts #
#----------------#
#these need to be the length of the science_frame list 
airmass_lis 	= [1.0,]
seeing_lis 		= [3.0,]

sky_sampling 	= 2.0 			#[float] Regridded sample size on sky in arcsec.
max_distance	= 102.0			#[float] Samples further away will have weight=0 [arcsec]. 
cube_sigma		= 4.8			#[float] Gaussian sigma of interpolation kernel [arcsec].
diffAtmRef		= True 			#[True/False] Differential atmospheric refraction correction applied if True.

#--------------------#
# collapse cube opts #
#--------------------#
col_wave_range 	= [0,0]			#[floatarray] Choose the wavelength range (IN ANGSTROMS) that you would like to build and image of 
									#you can choose [0,0] if you would like the collapse the entire data cube 


