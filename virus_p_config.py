############################################
# Configuration file for VIRUS-P reduction # 
############################################

#reduction_wrapper_virusP.py should be run inside of a folder containing the date folder containing raw data

redux_dir        = "VP_apr15_B1"	 	#name of the folder reduction is run 
date_folder      = "20150417"		#date folder containing raw data 

#Choose the spectral range of your VIRUS-P data 
spec_range   	 = [3500,5800] 		#Spectral range of that nights observations in angstroms

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
standard_frames  = ['Feige34','Feige67','Faige34']
#If you took sky frames, this should be a list of the object names of those frames as written in the headers
sky_frames		 = ['M82_sky','NGC4666_sky']
#This should be the object names of the science frames as written in the headers
science_frames      = ['M82_F1','M82_F2','NGC4666']	

#####################
# REDUCTION OPTIONS #
#####################
basic           = False		#run basic reduction (normalize, build mastertrace and masterarc)
run_deformer    = False		#run deformer to map spectral traces and build wavelength solution for fiber extraction
subsky          = False  	#run sky subtraction on sci images - Need to have run deformer, only runs on non-extracted spectra
fiberextract    = True	 	#extract spectra and save into fits file - Need to have run deformer
makecube 		= False		#builds data cube out of fiber extracted image - Need to have run fiberextract

################
# defomer opts #
################
ref_line		 = 5 			#[integer] arc line in line list that is the brightest but not blended (used as reference line)

########################
# sky subtraction opts #
########################
window_size 	= 250 			#[integer] Size (in image lines) of moving window for sky median
sky_kappa 		=[3.5,3.5]		#[floatarray] Lower and upper kappa for final sky kappa-sigma clipping.
smoothing 		= 5.0 			#[float] Smoothing factor for approximating spline.
sn_thresh		= 5 			#[float] Minimum signal to noise to flag fiber as continuum and ignore it during sky generation.

######################
# fiber extract opts #
######################
wl_resample 	 = True 		#[True/False] If True it will resample in wavelength, Does not resample in wavelength if False 
use_ap_corr	 	 = False		# Apply aperture correction using the fiber model (slow).

##################
# make cube opts #
##################
sky_sampling 	= 0.3 			#[float] Regridded sample size on sky in arcsec.
max_distance	= 5.0 			#[float] Samples further away will haver weight=0 [arcsec]. 
cube_sigma		= 0.75			#[float] Gaussian sigma of interpolation kernel [arcsec].
diffAtmRef		= True 			#[True/False] Differential atmospheric refraction correction applied if True.
