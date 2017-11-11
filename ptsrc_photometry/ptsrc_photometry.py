#! /usr/bin/env python

"""Performs point-source photometry using the ``IRAF`` tasks ``DAOFind``
and ``Phot`` on a single WFC3/UVIS image or list of images.
    
Used in the WFC3/UVIS contamination and stability monitoring program.    

Authors:

	C. Shanahan, Jul. 2016

    C.M. Gosmeyer, Nov. 2013

    D.M. Hammer, 2012
    
Use:

    Set appropriate parameters in the 'set_paths.py' script and in '__main__'

    >>> ssbx
    >>> python ptsrc_photom_uvis.py
    
Outputs:

    Ascii file. A photometry catalog named ``<filtername>_photcat.dat``.

History:

    Sept. 2012: 
        * Original script (v0.1).
    Oct. 2012: 
        * Added more robust handling of images with single chip 
          (e.g., IR, calibration).
    Nov. 2013: 
    	* Taken over by C.M. Gosmeyer (v0.3)
    Jan. 2014: 
        * Edited so that .coo and .mag files are saved.
        * Added function to retrieve name of filter.
        * Alphabetized the functions.
        * Deleted (now out-versioned) :func:`get_uvis_zeropoint` function.
        * Replaced `ptsrc_photom_uvis.py` version of 
          :func:`make_PAMcorr_image` with the version in 
          `ptsrc_photom_ircontam.py`.
    Feb 2014: 
        * DS9 now pops up if DAOfind finds more than one source and 
          allows user to select the source manually.
    Mar. 2014:
        * Deleted :func:`make_counts_image_uvis`, since 
          :func:`make_PAMcorr_image` seems to have replaced it and  
          just was never removed.  
    Jul. 2016:
    	*Taken over by Clare Shanahan
    	*Made some changes that make this script compatible with
    	 new file structure (including switch for temp directory)
    	*Implemented multiprocessing for speed
    	*Removed old and unused functionality to clean up code
    	*Deprecated the image select mode 
    	
    
Future Improvements:

    * Move from IRAF to photutils!  (see ``detector_tools/photometry``
      repo)
    * Store outputs in a database instead? (The structure for the 
      text file writing and reading is already there, clunky as it is; 
      may not be worth the time to refactor all the scripts even if 
      a database would be nice...)

Notes:

    * Be sure to change `pamdir` in :func:`make_PAMcorr_image` function to  
      your own directory containing the PAM correction FITS files.
    * If DS9 is having trouble, and you are in Ureka, it may be that the 
      Ureka DS9 is out-of-date.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')

import fileinput
import glob
import os
import pylab
import pyraf
import scipy
import shutil
from set_paths import *
from move_files_from_temp_directories import *

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from string import find
from scipy import interpolate
from time import sleep
from pyraf import iraf
from iraf import noao, digiphot, daophot
from astropy.table import Table
from astropy.io import fits as pyfits
from astropy.io import ascii
from stwcs.wcsutil import HSTWCS
from uvis_contam_tools import get_filter_name, get_wfc3_zeropoint, \
                              get_pivot_wavel, meanclip, make_PAMcorr_image
                              
from multiprocessing import Pool


#-------------------------------------------------------------------------------#

def append_phot_catalog(filt, perm_path_to_photocat,path_to_flts):
    """Appends new phot catalog to old catalog. Then need
    only run photometry over newest files, instead of all the 
    files every time.
    
    Parameters:
        filt : string
            Name of the filt
        origin : string
            Path to the photcat file. If '', assume in 
            current working directory.
            
    Returns:
        nothing
        
    Outputs:
        nothing
    """
    
    #if the photcat doesnt already exist, make a file and add the header
    if not os.path.isfile(perm_path_to_photocat+'/'+filt+"_photcat.dat"):
    	print "permenant photcat doesn't exist, making file for " + filt
    	with open(perm_path_to_photocat+'/'+filt+"_photcat.dat", "a") as old_photcat:
    		new_photcat = open(path_to_flts+'/'+filt+"_photcat.dat", "r")
    		#look in first few lines for # '
    		header_search = new_photcat.readlines()[0:3]
    		new_photcat.close()
    		for line in header_search:
    			if '#' in line:  
    				#print line  		
    				old_photcat.write(line)
    		
    		
    	
    	
    
    with open(perm_path_to_photocat+'/'+filt+"_photcat.dat", "a") as old_photcat:
        new_photcat = open(path_to_flts+'/'+filt+"_photcat.dat", "r")
        new_photcat_list = new_photcat.readlines()
        new_photcat.close()
        
        for line in new_photcat_list:
            if '#' not in line:
                old_photcat.write(line)


#-------------------------------------------------------------------------------#

def circular_mask(arr_shape, r, x_offset=0, y_offset=0):
    """Generates circular mask for 2D image.
        
    Parameters: 
        arr_shape : tuple of int
            Shape of the array to use the mask.
        r : int
            Radius of the mask in pixels.
        x_offset, y_offset : int or float, optional
            Mask offset relative to image center.
        
    Returns: 
        Numpy indices of the mask, rounded to nearest integer.
        
    Outputs:
        nothing
        
    References: 
        http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054470.html
    """
    assert len(arr_shape) == 2, 'Image is not 2-D'
    
    ny, nx = arr_shape
    assert nx > 1 and ny > 1, 'Image is too small'
    
    assert isinstance(r, (int, long)) and r > 0, 'Radius must be int > 0'
    
    xcen = np.round(0.5 * nx - 0.5 + x_offset).astype('int')
    ycen = np.round(0.5 * ny - 0.5 + y_offset).astype('int')
    
    x1, x2 = xcen - r, xcen + r
    y1, y2 = ycen - r, ycen + r
    
    assert y1 >= 0 and y2 < ny and x1 >= 0 and x2 < nx, \
          'Mask falls outside image bounds'
    
    y, x = np.ogrid[-r:r, -r:r]
    i = np.where(x**2 + y**2 <= r**2)
    
    a = np.zeros(arr_shape).astype('bool')
    a[y1:y2, x1:x2][i] = True
    
    return np.where(a)


#-------------------------------------------------------------------------------#

def get_modelpsf_uvis(wave_eval, rad_eval=[-9999.0]):
    """Retrieves model PSF. 
    
    Parameters:
        wave_eval : scalar of float
            Single wavelength.
        rad_eval : scalar of float
            Radius value.
            
    Returns:
        [rad_eval, flux_eval] : list
            List of radius value and flux value.  
    
    Outputs:
        nothing
    
    References:
        G. Hartig 2009 ISRs, http://www.stsci.edu/hst/wfc3/documents/ISRs/ 
    """
    # Checks
    if len(np.array([wave_eval])) != 1: 
        raise Exception('PSF may be evaluated at only 1 wavelength.')
    
    # read UVIS EE vs radius from Hartig
    extfile = 'Hartig_EE_model.dat'
    alldata = np.loadtxt(extfile)
    wave = alldata[0,1:]
    aper_rad = alldata[1:,0]
    
    # Match input wavelength to table wavelength (MUST MATCH)
    gdw = np.where(wave_eval == wave)[0]
    if len(gdw) == 1: data = alldata[1:,(gdw[0]+1)]
    else: raise Exception('Unique wavelength not found in external table.')
    
    # set radius positions to be evaluated & perform boundary checks
    nrads = data.size
    rmin=aper_rad[0]
    rmax=aper_rad[nrads-1]
    
    if rad_eval[0] < -999.0: #default is to use George's radius values
        rad_eval = aper_rad 
    if len(np.array(rad_eval)[rad_eval > rmax]) > 0 or \
        len(np.array(rad_eval)[rad_eval < rmin]) > 0: 
        raise Exception('Requested aperture radius is outside table boundaries.')
    
    # Evaluate table values at requested radii
    
    # Establish class for spline fitted model:
    interp = scipy.interpolate.InterpolatedUnivariateSpline(aper_rad, data)   
    if len(rad_eval) == 1: 
        flux_eval = np.array([interp(rad_eval)])
    else: 
        flux_eval = interp(rad_eval)
    
    # Check for wave input parameters outside exisiting boundaries or rules (OLD METHOD FOR ISR TABLE)
    #ny,nx=data.shape
    #xmin=wave[0]
    #xmax=wave[nx-1]
    #xeval = [wave_eval for x in xrange(len(rad_eval))]
    #if (wave_eval > xmax) or (wave_eval < xmin): raise Exception('Requested wavelength is outside table boundaries.')
    # Evaluate table at requested wave/radius by interpolating (APPLIES TO TABLE FROM ISR--NOT NEW TABLE)
    #interp = scipy.interpolate.RectBivariateSpline(aper_rad, wave, data)        # establish class for spline fitted model
    #flux_eval = interp.ev(rad_eval, xeval)
    
    return [rad_eval, flux_eval]


#-------------------------------------------------------------------------------#

def make_photom_catalog_uvis(data, filt, origin=''):
	"""Makes the photometry catalog ascii text file for UVIS data.
	
	Parameters:
	    data : dictionary
	        The photometry data for source star.
	    filt : string
	        Name of the filter.
        origin : string
            Path to the photcat file. If '', assume in 
            current working directory.
	        
	Returns:
	    nothing
	        
	Outputs:
	     ascii file. Photometry catalog named ``<filt>_photcat.dat``.
    """
	fnarr = [data[key]['filename'] for key in data.keys()]
	# Make sure filename does not include path.
	if '/' in fnarr[0]:
	    if 'temp_lacos' in fnarr[0]:
	        for i in range(len(fnarr)):
	            file_name = fnarr[i].split('temp_lacos/')[1]
	            fnarr[i] = file_name
	    else:
	        for i in range(len(fnarr)):
	            file_name = (fnarr[i].split('/'))[len(fnarr[i].split('/'))-1]
	            fnarr[i] = file_name
	
	amparr = [data[key]['amp'] for key in data.keys()]
	shutarr = [data[key]['shutter'] for key in data.keys()]
	mjdarr = [data[key]['mjd_avg'] for key in data.keys()]
	mjddeltarr = [data[key]['mjd_deltat'] for key in data.keys()]
	chiparr = [data[key]['chip'] for key in data.keys()]
	axis1arr = [data[key]['axis1'] for key in data.keys()]
	axis2arr = [data[key]['axis2'] for key in data.keys()]
	xcarr = [data[key]['xc'] for key in data.keys()]
	ycarr = [data[key]['yc'] for key in data.keys()]
	xcparr = [data[key]['xcp'] for key in data.keys()]
	ycparr = [data[key]['ycp'] for key in data.keys()]
	backarr = [data[key]['background'] for key in data.keys()]
	backrmsarr = [data[key]['background_rms'] for key in data.keys()]
	exptimearr = [data[key]['exptime'] for key in data.keys()]
	f1 = [data[key]['flux'][0] for key in data.keys()]
	f2 = [data[key]['flux'][1] for key in data.keys()]
	f3 = [data[key]['flux'][2] for key in data.keys()]
	f4 = [data[key]['flux'][3] for key in data.keys()]
	f5 = [data[key]['flux'][4] for key in data.keys()]
	f6 = [data[key]['flux'][5] for key in data.keys()]
	f7 = [data[key]['flux'][6] for key in data.keys()]
	f8 = [data[key]['flux'][7] for key in data.keys()]
	f9 = [data[key]['flux'][8] for key in data.keys()]
	f10 = [data[key]['flux'][9] for key in data.keys()]
	f12 = [data[key]['flux'][10] for key in data.keys()]
	f14 = [data[key]['flux'][11] for key in data.keys()]
	f16 = [data[key]['flux'][12] for key in data.keys()]
	f18 = [data[key]['flux'][13] for key in data.keys()]
	f20 = [data[key]['flux'][14] for key in data.keys()]
	f24 = [data[key]['flux'][15] for key in data.keys()]
	f28 = [data[key]['flux'][16] for key in data.keys()]
	f32 = [data[key]['flux'][17] for key in data.keys()]
	f36 = [data[key]['flux'][18] for key in data.keys()]
	f40 = [data[key]['flux'][19] for key in data.keys()]
	f45 = [data[key]['flux'][20] for key in data.keys()]
	f50 = [data[key]['flux'][21] for key in data.keys()]
	f55 = [data[key]['flux'][22] for key in data.keys()]
	f60 = [data[key]['flux'][23] for key in data.keys()]
	f65 = [data[key]['flux'][24] for key in data.keys()]
	f70 = [data[key]['flux'][25] for key in data.keys()]
    
	m1 = [data[key]['mag'][0] for key in data.keys()]
	m2 = [data[key]['mag'][1] for key in data.keys()]
	m3 = [data[key]['mag'][2] for key in data.keys()]
	m4 = [data[key]['mag'][3] for key in data.keys()]
	m5 = [data[key]['mag'][4] for key in data.keys()]
	m6 = [data[key]['mag'][5] for key in data.keys()]
	m7 = [data[key]['mag'][6] for key in data.keys()]
	m8 = [data[key]['mag'][7] for key in data.keys()]
	m9 = [data[key]['mag'][8] for key in data.keys()]
	m10 = [data[key]['mag'][9] for key in data.keys()]
	m12 = [data[key]['mag'][10] for key in data.keys()]
	m14 = [data[key]['mag'][11] for key in data.keys()]
	m16 = [data[key]['mag'][12] for key in data.keys()]
	m18 = [data[key]['mag'][13] for key in data.keys()]
	m20 = [data[key]['mag'][14] for key in data.keys()]
	m24 = [data[key]['mag'][15] for key in data.keys()]
	m28 = [data[key]['mag'][16] for key in data.keys()]
	m32 = [data[key]['mag'][17] for key in data.keys()]
	m36 = [data[key]['mag'][18] for key in data.keys()]
	m40 = [data[key]['mag'][19] for key in data.keys()]
	m45 = [data[key]['mag'][20] for key in data.keys()]
	m50 = [data[key]['mag'][21] for key in data.keys()]
	m55 = [data[key]['mag'][22] for key in data.keys()]
	m60 = [data[key]['mag'][23] for key in data.keys()]
	m65 = [data[key]['mag'][24] for key in data.keys()]
	m70 = [data[key]['mag'][25] for key in data.keys()]
	
	m1_err = [data[key]['merr'][0] for key in data.keys()]
	m2_err = [data[key]['merr'][1] for key in data.keys()]
	m3_err = [data[key]['merr'][2] for key in data.keys()]
	m4_err = [data[key]['merr'][3] for key in data.keys()]
	m5_err = [data[key]['merr'][4] for key in data.keys()]
	m6_err = [data[key]['merr'][5] for key in data.keys()]
	m7_err = [data[key]['merr'][6] for key in data.keys()]
	m8_err = [data[key]['merr'][7] for key in data.keys()]
	m9_err = [data[key]['merr'][8] for key in data.keys()]
	m10_err = [data[key]['merr'][9] for key in data.keys()]
	m12_err = [data[key]['merr'][10] for key in data.keys()]
	m14_err = [data[key]['merr'][11] for key in data.keys()]
	m16_err = [data[key]['merr'][12] for key in data.keys()]
	m18_err = [data[key]['merr'][13] for key in data.keys()]
	m20_err = [data[key]['merr'][14] for key in data.keys()]
	m24_err = [data[key]['merr'][15] for key in data.keys()]
	m28_err = [data[key]['merr'][16] for key in data.keys()]
	m32_err = [data[key]['merr'][17] for key in data.keys()]
	m36_err = [data[key]['merr'][18] for key in data.keys()]
	m40_err = [data[key]['merr'][19] for key in data.keys()]
	m45_err = [data[key]['merr'][20] for key in data.keys()]
	m50_err = [data[key]['merr'][21] for key in data.keys()]
	m55_err = [data[key]['merr'][22] for key in data.keys()]
	m60_err = [data[key]['merr'][23] for key in data.keys()]
	m65_err = [data[key]['merr'][24] for key in data.keys()]
	m70_err = [data[key]['merr'][25] for key in data.keys()]
    
	tt = {'#filename':fnarr, 'amp':amparr, 'shutter':shutarr, \
	      'mjd_avg':mjdarr, 'mjd_deltat':mjddeltarr, 'chip':chiparr, \
	      'axis1':axis1arr, 'axis2':axis2arr,'xc':xcarr, 'yc':ycarr, \
	      'xcp':xcparr, 'ycp':ycparr, 'background':backarr, \
	      'background_rms':backrmsarr, 'exptime':exptimearr, \
          'f1':f1, 'f2':f2, 'f3':f3,'f4':f4,'f5':f5,'f6':f6,'f7':f7,'f8':f8,\
          'f9':f9,'f10':f10,'f12':f12,'f14':f14,'f16':f16,'f18':f18,'f20':f20,\
          'f24':f24,'f28':f28,'f32':f32,'f36':f36,'f40':f40,'f45':f45,\
          'f50':f50,'f55':f55,'f60':f60,'f65':f65,'f70':f70,'m1':m1, 'm2':m2, \
          'm3':m3,'m4':m4,'m5':m5,'m6':m6,'m7':m7,'m8':m8,'m9':m9,'m10':m10,\
          'm12':m12,'m14':m14,'m16':m16,'m18':m18,'m20':m20,'m24':m24,\
          'm28':m28,'m32':m32,'m36':m36,'m40':m40,'m45':m45,'m50':m50,\
          'm55':m55,'m60':m60,'m65':m65,'m70':m70,'m1err':m1_err, \
          'm2err':m2_err, 'm3err':m3_err,'m4err':m4_err,'m5err':m5_err,\
          'm6err':m6_err,'m7err':m7_err,'m8err':m8_err,'m9err':m9_err,\
          'm10err':m10_err,'m12err':m12_err,'m14err':m14_err,'m16err':m16_err,\
          'm18err':m18_err,'m20err':m20_err,'m24err':m24_err,'m28err':m28_err,\
          'm32err':m32_err,'m36err':m36_err,'m40err':m40_err,'m45err':m45_err,\
          'm50err':m50_err,'m55err':m55_err,'m60err':m60_err,'m65err':m65_err,\
          'm70err':m70_err}

	ascii.write(tt, origin+filt+'_photcat.dat', \
	            names=['#filename','amp','shutter','mjd_avg','mjd_deltat',\
	                   'chip','axis1','axis2','xc','yc','xcp','ycp',\
	                   'background','background_rms','exptime', \
                       'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10',\
                       'f12','f14','f16','f18','f20','f24','f28','f32','f36',\
                       'f40','f45','f50','f55','f60','f65','f70',\
                       'm1','m2','m3','m4','m5','m6','m7','m8','m9','m10',\
                       'm12','m14','m16','m18','m20','m24','m28','m32','m36',\
                       'm40','m45','m50','m55','m60','m65','m70','m1err',\
                       'm2err','m3err','m4err','m5err','m6err','m7err',\
                       'm8err','m9err','m10err','m12err','m14err','m16err',\
                       'm18err','m20err','m24err','m28err','m32err','m36err',\
                       'm40err','m45err','m50err','m55err','m60err','m65err',\
                       'm70err'], \
                formats={'#filename':'%s','amp':'%s','shutter':'%s',\
                         'mjd_avg':'%9.4f','mjd_deltat':'%6.4f','chip':'%i',\
                         'axis1':'%i','axis2':'%i','xc':'%8.3f','yc':'%8.3f',\
                         'xcp':'%8.3f','ycp':'%8.3f', 'background':'%0.5f',\
                         'background_rms':'%0.5f', 'exptime':'%0.2f', \
                         'f1':'%0.3f', 'f2':'%0.3f','f3':'%0.3f','f4':'%0.3f',\
                         'f5':'%0.3f','f6':'%0.3f','f7':'%0.3f','f8':'%0.3f',\
                         'f9':'%0.3f','f10':'%0.3f','f12':'%0.3f',\
                         'f14':'%0.3f','f16':'%0.3f','f18':'%0.3f',\
                         'f20':'%0.3f','f24':'%0.3f','f28':'%0.3f',\
                         'f32':'%0.3f','f36':'%0.3f','f40':'%0.3f',\
                         'f45':'%0.3f','f50':'%0.3f','f55':'%0.3f',\
                         'f60':'%0.3f','f65':'%0.3f','f70':'%0.3f',\
                         'm1':'%0.3f','m2':'%0.3f','m3':'%0.3f','m4':'%0.3f',\
                         'm5':'%0.3f','m6':'%0.3f','m7':'%0.3f','m8':'%0.3f',\
                         'm9':'%0.3f','m10':'%0.3f','m12':'%0.3f',\
                         'm14':'%0.3f','m16':'%0.3f','m18':'%0.3f',\
                         'm20':'%0.3f','m24':'%0.3f','m28':'%0.3f',\
                         'm32':'%0.3f','m36':'%0.3f','m40':'%0.3f',\
                         'm45':'%0.3f','m50':'%0.3f','m55':'%0.3f',\
                         'm60':'%0.3f','m65':'%0.3f','m70':'%0.3f', \
                         'm1err':'%0.3f', 'm2err':'%0.3f','m3err':'%0.3f',\
                         'm4err':'%0.3f','m5err':'%0.3f','m6err':'%0.3f',\
                         'm7err':'%0.3f','m8err':'%0.3f','m9err':'%0.3f',\
                         'm10err':'%0.3f','m12err':'%0.3f','m14err':'%0.3f',\
                         'm16err':'%0.3f','m18err':'%0.3f','m20err':'%0.3f',\
                         'm24err':'%0.3f','m28err':'%0.3f','m32err':'%0.3f',\
                         'm36err':'%0.3f','m40err':'%0.3f','m45err':'%0.3f',\
                         'm50err':'%0.3f','m55err':'%0.3f','m60err':'%0.3f',\
                         'm65err':'%0.3f','m70err':'%0.3f'})


#-------------------------------------------------------------------------------#

def make_source_location_histogram_plots_uvis(data, file_name, ff, im, coordfile, \
    filt, path_to_cleans=''):
	"""Makes the diagnostic source location and histogram plots for 
	UVIS data.
	
	Parameters:
	    data : dictionary 
	        The photometry data for source star.
	    file_name : string
	        Name of the FITS file. 
	    ff : int
	        Index of the file in a list of the files in the directory.
	    im : array
	        The image array from the FITS file. 
        coordfile : string
            `coords` param for ``IRAF/Phot``. Will be ``0.coo.1`` file if 
            from :func:`run_daofind``.
	    filt : string
	        Name of the filter.
	    path_to_cleans : string
	        Path to the ``*clean.fits`` and associated files. If '', assumes
	        current working directory. 
	       
	Returns:
	    nothing       
	        
	Outputs:
	    PNG file. Plots of the object's location and the histogram of 
	    its background in file called ``<file rootname>_srcloc.png``. 
    """
	pylab.ion()
	if ff == 0:
		fig = pylab.figure()
		fig.subplots_adjust(wspace=0.4)
	else:
		pylab.clf()
		
	xc,yc = np.loadtxt(coordfile, unpack=True, usecols = (0,1)) 
	# plot #1 - object position
	sz=50.0
	x0=np.round(xc)-sz/2.
	x1=np.round(xc)+sz/2.
	y0=np.round(yc)-sz/2.
	y1=np.round(yc)+sz/2.
	ax1 = pylab.subplot(1,2,1)
	ax1.imshow(np.log10(im[y0:y1,x0:x1]),interpolation='nearest')
	ax1.autoscale(axis='both',enable=False)
	ax1.scatter([xc-x0-1.0], [yc-y0-1.0], marker='x', s=200., color='w')
	pylab.title('X = '+str(xc)+'  Y = '+str(yc))

	# plot #2 - background histogram
	tmp_image=glob.glob(path_to_cleans + '*back.fits')[0]
	backim = pyfits.getdata(tmp_image)
	#--measure back statistics (mean and mode via IRAF)
	initback = iraf.imstatistics(tmp_image+'[0]', fields='mode,stddev', \
	                             lower = -100, upper = 10000, nclip=7, \
	                             lsigma=3.0, usigma=3.0, cache='yes', \
	                             format='no',Stdout=1)
	#print 'initback:'
	#print initback
	if 'INDEF' not in initback[0]:
		llim = float(initback[0].split('  ')[0]) - 10.0*\
				float(initback[0].split('  ')[1])
		ulim = float(initback[0].split('  ')[0]) + 10.0*\
	           float(initback[0].split('  ')[1])
		backstats=iraf.imstatistics(tmp_image+'[0]', fields='mean,mode', \
	                               lower=llim, upper=ulim, nclip=7,lsigma=3.0, \
	                               usigma=3.0, cache='yes', format='no',Stdout=1)
		backmean=float(backstats[0].split('  ')[0])
		backmode=float(backstats[0].split('  ')[1])
		fbackim= np.ndarray.flatten(backim)
		gd=np.where((fbackim > llim) & (fbackim < ulim))[0]
		backmedian=meanclip(fbackim[gd],maxiter=7,return_median=1)[0]

		ax2 = pylab.subplot(1,2,2)
		pylab.hist(fbackim[gd],log=True)
		pylab.ylim(0.5,600000)
		pylab.xlim(-20,20)
		pylab.plot([backmode,backmode],[0.5,600000],ls='-',color='red',\
	               label='mode')
		pylab.plot([backmedian,backmedian],[0.5,600000],ls='--',color='aqua',\
    	           label='median')
		pylab.plot([backmean,backmean],[0.5,600000],ls=':',color='black',\
    	           label='mean')
		pylab.legend(loc=2, handletextpad=0.0, borderpad=0.0, frameon=False, \
    	             handlelength=1.)
		pylab.title('Histogram of Background Pixels')
		pylab.xlabel('Background [e-]')
		pylab.ylabel('Number of Objects')
		pylab.annotate('chip '+str(data[ff]['chip']), [0.77,0.95], \
    	               xycoords='axes fraction')
		pylab.annotate(filt,[0.77,0.80],xycoords='axes fraction')

		
	pylab.savefig(file_name.split('.fits')[0]+'_srcloc.png')
	pylab.ioff()


#-------------------------------------------------------------------------------#

def move_files(origin=''):
	"""Moves position-histogram PNG plots into subdirectory.
	
	These are the plots created by the function
	:func:`make_source_location_histogram_plots_uvis`.
	
	Parameters:
        origin : string
            Path to the photcat file. If '', assume in 
            current working directory.
	    
	Returns:
	    nothing
	    
	Outputs:
	    nothing
    """
	png_file_list = glob.glob(origin+'*png')
	if png_file_list != []:
		if not os.path.exists(origin+'positions-histograms'):
			os.makedirs(origin+'positions-histograms')
		for png in png_file_list:
			shutil.move(str(png), origin+'positions-histograms')

#-------------------------------------------------------------------------------#

def open_ds9_imexam(file_name, ext, coordfile):
    """Opens DS9 window and sets up ``IRAF/IMEXAM`` so that user can
    manually select location of star. The location is stored in a 
    numpy array.
    
    Parameters: 
        file_name : string
            Name of *_FLT.CLEAN.FITS.
        ext : integer
            Extension of image.
        coordfile : string
            Name of the *coo.1 file.
            
    Returns:
        
    Outputs:
        nothing
    """
    print "--> After select frame, center cursar on star and click."
    print "--> You are now in iraf.imexam. Type 'm' and then 'q'."
    ds9_window = os.spawnlp(os.P_NOWAIT,"ds9", "ds9", "-title", "ptsrc_photom")
    sleep(2) # DS9 needs nap, else grumpy
    iraf.display(file_name+ "[" + str(ext) + "]")
    iraf.imexam(logfile="imexam.txt", keeplog="yes")
    x_ds9, y_ds9 = readlogfile("imexam.txt")
    print "imexam x, y coords: "
    print x_ds9, y_ds9
    print "Closing DS9 window..."
    os.kill(ds9_window, 1)	# Seems to be working; I hope it is closing the window safely...
    np.savetxt(coordfile, zip(np.array([x_ds9]), \
               np.array([y_ds9])), fmt='%0.5f')
    os.remove("imexam.txt")
    
    #return find_name


#-------------------------------------------------------------------------------#

def readlogfile(file_name):
    """Reads in the reference star locations from the log file.
    
    Parameters:
        file_name : text file
            Contains the source's x and y positions, etc.
    
    Returns:
        x[0], y[0] : floats
            The source's x and y positions only.
    Outputs:
        nothing
    
    References: 
        from http://www.astronomy.ohio-state.edu/~martini/osmos/oalign.py
    """
    F = open(file_name, "r")
    M = F.readlines()
    F.close()
    x = []
    y = []
    
    for i in range(len(M)):
        if find(M[i], "#") <0:
            # Since the x,y are recorded as, e.g., '[398:402,300:304]'
            x.append(float(M[i].split('[')[1].split(':')[0]) + 2.0)
            y.append(float(M[i].split(',')[1].split(':')[0]) + 2.0)
    
    # Notice this function allows multiple stars to be chosen. For our purposes here,
    # we need pick only one, and so return the first values in the list.
    return x[0], y[0]


#-------------------------------------------------------------------------------#

def rename_photcat(filt, origin='', revert=True):
    """Renames the ``<filt>_photcat.dat`` file
    to or back from ``<filt>_photcat.store.dat``.
    
    Parameters:
        filt : string
            Name of the filter.
        origin : string
            Path to the photcat file. If '', assume in
            current working directory.
        revert : {True, False}
            If True, renames photcat.store.dat back to photcat.dat.
            If False, renames photcat.dat file to photcat.store.dat.
            
    Returns:
        nothing
        
    Outputs:
        nothing
    """
    if revert == False:
        os.rename(origin+filt+"_photcat.dat", origin+filt+\
                  "_photcat.store.dat")
        
    if revert == True:
        os.rename(origin+filt+"_photcat.store.dat", origin+filt+\
                  "_photcat.dat")
        

#-------------------------------------------------------------------------------#

def replace_filevalue(file_name, orgval, newval):
    """Replaces unwanted values in external file.
    
    Parameters:
        file_name : text file
            The input file, assumed to be a ``coo.1`` (coordinate) file.
        orgval : int, float, or string
            Original value.
        newval : int, float, or string
            New value.
            
    Returns:
        nothing
            
    Outputs:
        The modified input file. 
    """
    for line in fileinput.input(file_name, inplace = 1):
        print line.replace(str(orgval), str(newval)),
    fileinput.close()


#-------------------------------------------------------------------------------#

def run_daofind(image, extension=0, outfile='default',dthreshold=3.0, \
    fwhmpsf=2.5, backsigma=-1.0,rdnoise=-1.0):
    """Runs ``IRAF/DAOFind`` on input image.
    
    Parameters:
        image : FITS file 
            Input image, assumed to be a WFC3/UVIS or WFC3/IR FITS file.
            Can handle UVIS subarrays. 
        extension : int
            Use `extension = 1` ('SCI') for unaltered FLT. Use 
            `extension = 0` for cosmic-ray cleaned FLT.
        outfile : string
            By default, the outfile is named ``<image>0.coo.1``.
        dthreshold : float
            `threshold` param for ``IRAF/DAOFind``.
        fwhmpsf : float
            `fwhmpsf` param for ``IRAF/DAOFind``.
        backsigma : float
            `sigma` param for ``IRAF/DAOFind``.
        rdnoise : float
            Used to calculate `readnoise` param for ``IRAF/DAOFind``.
            
    Returns:
        outfile : text file
            File containing the coordinates of sources found with  
            ``DAOFind``, named ``<image>_0.coo.1``.
            
    Outputs:
        nothing
    """
    # Parse input parameters
    if outfile == 'default':
        outfile = image+'0.coo.1'
    
    # Read in fits header
    f = pyfits.open(image)
    fheader = f[0].header
    
    # Extract relevant info from the header 
    # (exposure, filter, input/output pixel scale)
    exptime = fheader['exptime']
    instr = fheader['INSTRUME']
    if instr == 'WFC3':
        filt = fheader['FILTER']
    else: #assuming ACS
        filt = fheader['FILTER1']
        if filt[0] == 'C':
            filt == fheader['FILTER2']
    opxscl=0.03962
    ipxscl=0.03962
    f.close()

    # Assign number of flt images (IR/calibration images only have 1 chip; 
    # NDRIZIM keyword includes both chips from single FLT)
    if (fheader['detector'] == 'IR'):
        nchips = 1.0 		# IR
    elif (fheader['subarray'] == True) and (len(fheader['CCDAMP']) == 1):
    	nchips = 1.0        # UVIS sub-array
    elif (fheader['detector'] == 'UVIS') and (fheader['subarray'] == False):
        nchips = 2.0        # UVIS full-frame
    else:
        raise exception('Image type is not defined.')
    num_flts = 1.0

    # Perform read noise correction
    if rdnoise < 0.0:
        amps = fheader['CCDAMP']
        rdnoise = np.zeros(len(amps))
        for namp in xrange(len(amps)): 
            rdnoise[namp] = fheader['READNSE'+amps[namp]]
    rdnoise_corr = np.sqrt(num_flts * (np.average(rdnoise) * opxscl/ipxscl)**2)

    # Perform background noise calculation
    if backsigma < 0.0:
        backstats=iraf.imstatistics(image+'[0]', fields='stddev', \
                                    lower = -100, upper = 100, nclip=5, \
                                    lsigma=3.0, usigma=3.0, cache='yes', \
                                    format='no',Stdout=1)
        if backstats[0] == 'INDEF':
            backstats[0] = '0' 	#Added 13 Jan 2014
        backsigma=float(backstats[0])

    # remove old daofind files
    file_query = os.access(outfile, os.R_OK)
    if file_query == True:
        os.remove(outfile)
    iraf.daofind.unlearn()
    iraf.daofind(image=image+'[' + str(extension) + ']', interactive='no', \
            verify='no',output=outfile, fwhmpsf=fwhmpsf, \
            sigma=backsigma, readnoise=rdnoise_corr, itime=exptime, \
            threshold=dthreshold, datamin=-10, datamax=100000)

    # Display results of daofind (***WORK IN PROGRESS***)
    #os.system('ds9&')
    #tmp=image.split('_cnts')
    #iraf.display(tmp[0]+'.fits',1, zscale='no', zrange='no', z1=0, z2=100,ztrans='log')
    #iraf.tvmark(1,outfile,mark = 'circle', radii = 8, color = 205)

    return outfile                # return name of coordinate file


#-------------------------------------------------------------------------------#

def run_daophot_uvis(image, outfile='default', coordfile='NA', \
                     backmethod='mean', backval=None, backsigma=None, \
                     rdnoise=None, \
                     apertures='1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,45,50,55,60,65,70', \
                     cbox=3.0, annulus=17.0, dannulus=3.0, \
                     calgorithm='centroid', salgorithm='median', fwhmpsf=2.5, \
                     extension=0):
    """Runs ``IRAF/Phot`` on input UVIS image.
    
    Parameters: 
        image : string
            Name of the FITS image to be run photometry on.
        outfile : string
            `outfile` param for ``IRAF/Phot``. Named ``<image>0.mag`` by 
            default.
        coordfile : string
            `coords` param for ``IRAF/Phot``. Will be ``0.coo.1`` file if 
            from :func:`run_daofind``.
        backmethod : string
            The `mean`, `median`, and `mode` are all supported.
        backval : float
            If set to None, will be the background value found from the 
            `backmethod`.
        backsigma : float
            If set to None, will be the background rms found from the 
            `backmethod`.
        rdnoise : float
            Used to calculate `readnoise` param for ``IRAF/Phot``.
        apertures : string of number list
            `apertures` param for ``IRAF/Phot``. The aperture radii 
            (pixels) are ints separated by commas. 
        cbox : float
            `cbox` param for ``IRAF/Phot``.
        annulus : float
            param currently not in use
        dannulus : float
            param currently not in use
        calgorithm : string
            `calgorithm` param for ``IRAF/Phot``.
        salgorithm : string
            `calgorithm` param for ``IRAF/Phot``.
        fwhmpsf : float
            `fwhmpsf` param for ``IRAF/Phot``.
        extension : int
            Use `extension = 1` ('SCI') for unaltered FLT. Use 
            `extension = 0` for cosmic-ray cleaned FLT.
            
    Returns:
        backval, backsigma : floats
            The background (sky) stats for image.
    
    Outputs: 
        Text file. Contains output from ``Phot``, including source's
        fluxes for the specified aperture radii. Named ``<image>_0.mag.1``
        
    Example: 
        >>> back, backrms = run_daophot_uvis(cnts_name[ff], coordfile=find_name[ff], 
        >>>                                  outfile=phot_name[ff], calgorithm='gauss', 
        >>>                                  cbox=10., backmethod=backmeth, extension=0)
    """
    # Parse input parameters.
    if outfile == 'default':
        outfile = image + '0.mag'
    if coordfile == 'NA':
        coordfile = image + '0.coo'

    # Extract header info.
    prihdr = pyfits.getheader(image)
    exptime = prihdr['exptime']
    instrum = prihdr['INSTRUME']
    detector = prihdr['DETECTOR']
    SUBARRAY = prihdr['SUBARRAY']
    ccdamp = prihdr['CCDAMP']

    if instrum == 'WFC3':
        filt = prihdr['FILTER']
    elif instrum == 'ACS':
        filt = prihdr['FILTER1']
        if filt[0] == 'C':
            filt == prihdr['FILTER2']
    else:
        raise Exception('Instrument '+instrum+' not covered in our case list.')

    # Record native pixel scale and no. of chips.
    if instrum == 'WFC3':
        if detector == 'UVIS':
            pscale_nat = 0.03962
            if ((SUBARRAY == True) & (len(ccdamp) == 1)): nchips = 1.0
            elif SUBARRAY == False: nchips = 2.0
            else: raise Exception('Image type is not defined.')
        elif detector == 'IR':
            pscale_nat = 0.12825
            nchips = 1.0
        else:
            raise Exception('WFC3 Detector '+detector+' not covered in our case list.')
    elif instrum == 'ACS':
        if detector == 'WFC':
            pscale_nat = 0.049
            if ((SUBARRAY == True) & (len(ccdamp) == 1)):
                nchips = 1.0
            elif SUBARRAY == False:
                nchips = 2.0
            else:
                raise Exception('Image type is not defined.')
        else:
            raise Exception('ACS Detector '+detector+' not covered in our case list.')
    else:
        raise Exception('Instrument '+instr+' not covered in our case list.')


    # Record pixel scale of current image, image type, image axis lengths, 
    # and number of flts.
    sciext = []
    pscale_img = prihdr.get('D001SCAL',default='NA')
    if pscale_img == 'NA':
        imtype = 'flt'
        # we dont distinguish between flt/crclean, i.e., assume pscales are equal            
        pscale_img = pscale_nat
        num_flts = 1.0
        naxis1 = pyfits.getval(image,'NAXIS1',ext=extension)  #ext = ('SCI',1), for cleaned, ext = 0
        naxis2 = pyfits.getval(image,'NAXIS2',ext=extension)
        # -- record location of science extension
        hdulist = pyfits.open(image)
        for ext in xrange(len(hdulist)):
            if hdulist[ext].name == 'SCI': sciext.append(ext)
        hdulist.close()
        if len(sciext) != 1:
            sciext.append(0) #raise Exception('We do not handle images with '+str(len(sciext))+' SCI extensions.')
    else:
        imtype ='drz'
        num_flts = prihdr['NDRIZIM']/nchips
        naxis1 = pyfits.getval(image, 'NAXIS1', ext=extension)   #ext = ('SCI',1), for cleaned, ext = 0
        naxis2 = pyfits.getval(image, 'NAXIS2', ext=extension)
        sciext.append(0)

    # Get zeropoints.
    if instrum == 'WFC3':
        zeropt = get_wfc3_zeropoint(filt)
    elif instrum == 'ACS' and imtype == 'drz':
        zeropt = get_acs_zeropoint(prihdr)
    elif instrum == 'ACS' and imtype == 'flt':
        zeropt = get_acs_zeropoint(pyfits.getheader(image,ext=('SCI',1)))

    # Estimate read noise.
    if rdnoise == None:
        rdnoise = np.zeros(len(ccdamp))
        for namp in xrange(len(ccdamp)):
            rdnoise[namp] = prihdr['READNSE'+ccdamp[namp]]
    rdnoise_corr = np.sqrt(num_flts * (np.average(rdnoise) * \
                           pscale_img/pscale_nat)**2)


    # Measure the background and noise
    if ((backval == None) | (backsigma == None)):
        # Read in the x/y center of the source.
        xc,yc = np.loadtxt(coordfile, unpack=True, usecols = (0,1))
        
        # Create temporary image for bckgrd measurement that masks sources 
        # out to 80 pixels (assign a very low number).
        tmp_image = image+'.back.fits'
        shutil.copy(image, tmp_image)
        hdulist = pyfits.open(tmp_image, mode='update')
        maskim = hdulist[sciext[0]].data
        if detector == 'IR':
            maskrad = 30
        else:
            maskrad = 80
        maskim[circular_mask(maskim.shape, maskrad, x_offset=(xc-naxis1/2.0), \
               y_offset=(yc-naxis2/2.0))] = -99999.0

        # Also mask out sources with zero effective exposure 
        # [WE ELIMINATE PIXELS WITHIN 20 OF IMAGE BORDER].
        maskim[:,0:20] = -99999.0
        maskim[:,-20:] = -99999.0
        maskim[0:20,:] = -99999.0
        maskim[-20:,:] = -99999.0

        # Generate initial guess for lower/upper limits (use 10 sigma).
        fmaskim = np.ndarray.flatten(maskim)
        llim = -100
        ulim = 10000.0
        init_median,init_rms = meanclip(fmaskim[(fmaskim > llim) & \
                                        (fmaskim < ulim)], maxiter=7, \
                                        return_median=1)
        llim = init_median - 10.0*init_rms
        ulim = init_median + 10.0*init_rms

        # Measure background and rms.

        if backmethod.lower() == 'mean':
            back,backrms=meanclip(fmaskim[(fmaskim > llim) & \
                                  (fmaskim < ulim)], maxiter=7)
        elif backmethod.lower() == 'median':
            back,backrms = meanclip(fmaskim[(fmaskim > llim) & \
                                    (fmaskim < ulim)],maxiter=7,return_median=1)
        elif backmethod.lower() == 'mode':
            backmean,backrms = meanclip(fmaskim[(fmaskim > llim) & \
                                        (fmaskim < ulim)], maxiter=7)
            nbins = np.ceil(80.0/(0.1*backrms))
            cc,bb,pp = pylab.hist(fmaskim[(fmaskim > llim) & \
                                  (fmaskim < ulim)], log=True, bins=nbins, \
                                  range=(-40.0,40.0))
            back = bb[cc.argmax()] + (bb.max()-bb.min())/(2.0*(len(bb)-1))
        else:
            raise Exception('Background statistical method ' + backmethod + \
                            ' is not covered in our case list.')

        if backval == None:
            backval = back
        if backsigma == None:
            backsigma = backrms

        #print '\n BACKGROUND =  '+str(backval)
        #print ' BACKGROUND RMS =  '+str(backsigma)+' \n'


        # Case of no aperture size given (we select aperture sizes of 
        # UVIS=0.2/0.4", IR=0.27/0.4", ACS=0.25/0.5").
        if apertures == '':
            if instrum == 'WFC3' and detector == 'IR':
                apertures=str(0.27/pscale_img)+','+str(0.4/pscale_img)
            elif instrum == 'WFC3' and detector == 'UVIS':
                apertures=str(0.2/pscale_img)+','+str(0.4/pscale_img)
            elif instrum == 'ACS' and detector == 'WFC':
                apertures=str(0.25/pscale_img)+','+str(0.5/pscale_img)
            else:
                raise exception('Instrument/Detector '+ instrum + '/' \
                                 + detector + ' not covered in case list.')

        # Remove old phot output files.
        file_query = os.access(outfile, os.R_OK)
        if file_query == True:
            os.remove(outfile)

        # Run phot.
        iraf.phot.unlearn()         # reset daophot parameters to default values
        iraf.phot(image=image+'['+str(sciext[0])+']', interactive='no', \
                verify='no', coords=coordfile, output=outfile, \
                fwhmpsf=fwhmpsf, sigma=backsigma, readnoise=rdnoise_corr, \
                itime=exptime, calgorithm=calgorithm, cbox=cbox, \
                skyvalue=backval, apertures=apertures, zmag=zeropt, \
                salgorithm='constant')
        #annulus=annulus, dannulus=dannulus


        return backval, backsigma    # return computed background stats for image


#-------------------------------------------------------------------------------#
# The main controller, and then some.
# Creates the photometry catalog.
# Also creates sanity-check plots for each target star (each flt/flc file).
# (Might someday try to slice this up into more readable functions.)
#-------------------------------------------------------------------------------#

def ptsrc_photom_uvis(args):
    PAMcorr_switch=True
    images_seperate=False
    use_existing_coor=False
    force_ds9_select=False
    backmethod='mean'
    """Main for creating UVIS photometry catalog. 
    
    Runs source finding and photometry on all FLTs in given 
    directory `path_to_photcat`.
    
    Parameters:
        PAMcorr_switch : {True, False}
            On by default. Switch off if do NOT want to pixel area map
            correction performed on images.
        images_seperate : {True, False}
            Off by default. Switch on if the cleaned FLTs are in a 
            subdirectory `temp_lacos` instead of on same level as the 
            `<filt>_photcat.dat` file.
        use_existing_coor : {True, False}
            Off by default. Switch on if want to use any .coo files
            already in existence.
        force_ds9_select : {True, False}
            Off by default. Switch on if you want to force ds9 to open
            on every star so you can manually select them. 
        path_to_photcat : string
            Path to the photcat file. If '', assume in
            current working directory.
        backmethod : string
            The `mean`, `median`, and `mode` are all supported. Default of 'mean'.
    
    Returns:
        nothing
        
    Outputs:
        ascii file. `<filt>_photcat.dat`.
        0.coo.1 files. Contain source location (from ``DAOFIND``) for 
            each image.
        0.mag.1 files. Contain photometry (from ``PHOT``) for each 
            image.
        png files. `<rootname>_fl*.clean_srcloc.png`.
            Sanity check plot for each image: shows the location of 
            source found by ``PHOT`` and a 'brightness profile'. 
    """ 

    #unpack args  
    path_to_flts=args[0]
    path_to_photcat=path_to_flts #unfortunate naming...path_to_photcat is actually where flts are
    perm_path_to_photcat=args[1]
    
    append_to_photcat=True
    
    origin=path_to_photcat
    fits_list=glob.glob(path_to_flts+'/'+'*clean.fits')
    fits_list=fits_list
    if len(fits_list)==0:
    	print 'no flts found in ' + path_to_flts
    	return
    filt = pyfits.open(fits_list[0])[0].header['FILTER']
    
    
    # If a photcat file already exists, rename it so it will be saved.
    # After photometry is run, append the new photcat file to the saved
    # photcat file, and rename back to <filt>_photcat.dat.
    photcat_list = glob.glob(path_to_photcat + '*photcat.dat')
    if photcat_list != []:
        rename_photcat(filt, origin=path_to_photcat, revert=False)
    
    # Initialize filename and aperture correction variables
    if not images_seperate: 
        file_list = fits_list
    elif images_seperate:
        file_list = fits_list
    cnts_name = []
    for file_name, ff in zip(file_list, xrange(len(file_list))):
        tmp=file_name.split('.fits')
        cnts_name.append(tmp[0] + '_cnts.fits')
    
    find_name = [cnts_name[x]+'0.coo.1' for x in xrange(len(cnts_name))]
    phot_name = [cnts_name[x]+'0.mag.1' for x in xrange(len(cnts_name))]
    
    # Initialize a dictionary (structure) to hold data for each image
    data = {}
    
    # Generate source catalogs: src detection & photometry within same 
    # image (i.e., not color catalogs)
    for file_name, ff in zip(file_list,xrange(len(file_list))):
    
        print "--------------------------"
        print file_name
        print ff
            
        # Since LACosmic cleaned files have only 0th extension. 
        # Otherwise use SCI (1st) extension.
        if 'clean' not in file_name and '_dr' not in file_name:
            ext = 1
        else:
            ext = 0
        
        # Check whether to do PAM correction.
        if PAMcorr_switch:
            make_PAMcorr_image(file_name,outfile=cnts_name[ff], extension=0)
        
        # Check whether user wants to use already existing .coo file.
        if use_existing_coor and os.path.isfile(find_name[ff]): 
            print "USING PRE-EXISTING .COO FILE."
        
        # Or whether to force manual selection.
        elif force_ds9_select:
            open_ds9_imexam(file_name, ext, find_name[ff])
            replace_filevalue(find_name[ff], 'INDEF', 0.0)
        
        # Or whether to let DAOFIND do the job.
        else:

            if PAMcorr_switch:
                run_daofind(cnts_name[ff], outfile=find_name[ff], \
                            dthreshold=800.0, extension=ext)
        
            elif not PAMcorr_switch:
                run_daofind(file_name, outfile=find_name[ff], dthreshold=800.0, \
                            extension=ext)        
        
            # This try-except-else loop is for finding the source automatically, 
            # then manually with DS9 if all else fails.
            try:
                replace_filevalue(find_name[ff], 'INDEF', 0.0)
                xx,yy,mm,sharp,round,ground,id = np.loadtxt(find_name[ff], \
                                                            unpack=True)
            except ValueError:
                #print "Oops! Too little starlight found. Forcing DS9 selection..."
                #open_ds9_imexam(file_name, ext, find_name[ff])
                #replace_filevalue(find_name[ff], 'INDEF', 0.0)
                f = open(data_dir+'/bad_images.txt','a')  
                f.write(str(file_name)+'\n')
                f.close()
                print 'got a bad one'
                continue
            else:
                #Whatever. Try except won't work without an else?

				# Run through a decision tree to select the correct object if
				#  more than detected by find.
				if xx.size == 0:
					print 'NO OBJECTS WERE DETECTED.'
				elif xx.size > 1:
					# If more than 1 object detected - first remove any objects 
					# near borders as these are not expected
					xxx = xx[(xx > 150.) & (xx < 356.) & (yy > 150.) & (yy < 356.)]
					yyy = yy[(xx > 150.) & (xx < 356.) & (yy > 150.) & (yy < 356.)]
					mmm = mm[(xx > 150.) & (xx < 356.) & (yy > 150.) & (yy < 356.)]
					sxx = sharp[(xx > 150.) & (xx < 356.) & (yy > 150.) & (yy < 356.)]
					rxx = round[(xx > 150.) & (xx < 356.) & (yy > 150.) & (yy < 356.)]
					gxx = ground[(xx > 150.) & (xx < 356.) & (yy > 150.) & (yy < 356.)]
					if xxx.size == 0:
						print 'NO OBJECTS MATCHED IMAGE LOCATION CRITERIA.'
						# If there is more than one source:
						# Open DS9, manually select the star. The coords should save to the .coo file.
						open_ds9_imexam(file_name, ext, find_name[ff])
					elif xxx.size == 1:
						np.savetxt(find_name[ff],zip(xxx,yyy),fmt='%0.5f')
					else:
						# select the objects by find geometric parameters
						gd = np.where((sxx < 0.82) & (np.abs(gxx) < 0.4))[0]
						if gd.size == 1:
							np.savetxt(find_name[ff],zip(xxx[gd],yyy[gd]),fmt='%0.5f')
							sxx = sxx[gd]
							rxx = rxx[gd]
							gxx = gxx[gd]
						elif gd.size > 1:
							# select the object with the brightest magnitude
	#                        #dist = np.sqrt((xxx[gd]-255.0)**2 + (yyy[gd]-255.0)**2)
							mag = mmm[gd]
							xxx = xxx[gd][mag == np.min(mag)]
							yyy = yyy[gd][mag == np.min(mag)]
							sxx = sxx[gd][mag == np.min(mag)]
							rxx = rxx[gd][mag == np.min(mag)]
							gxx = gxx[gd][mag == np.min(mag)]
							np.savetxt(find_name[ff],zip(xxx,yyy),fmt='%0.5f')
						else:
							print 'NO OBJECTS MATCHED SHARPNESS/ROUNDNESS CRITERIA.'
							# If there is more than one source:
							# Open DS9, manually select the star. 
							# The coords should save to the .coo file.
							open_ds9_imexam(file_name, ext, find_name[ff])

				else:
					sxx = np.array([sharp])
					rxx = np.array([round])
					gxx = np.array([ground])
					


        # Once have found our source, do photometry:
        if PAMcorr_switch:
            back, backrms = run_daophot_uvis(cnts_name[ff], \
                            coordfile=find_name[ff], \
                            outfile=phot_name[ff], calgorithm='gauss', \
                            cbox=10., backmethod=backmethod, extension=ext)
        if not PAMcorr_switch:
            back, backrms = run_daophot_uvis(file_name, \
                                             coordfile=find_name[ff], \
                                             outfile=phot_name[ff], \
                                             calgorithm='gauss', cbox=10., \
                                             backmethod=backmethod, extension=ext)
        replace_filevalue(phot_name[ff], 'INDEF',-9999.0)
        iraf.txdump(phot_name[ff],'xcenter,ycenter,flux, mag, merr', \
                    'yes', Stdout=phot_name[ff]+'.trimmed')
        xc,yc,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,\
        f18,f19,f20,f21,f22,f23,f24,f25,f26,m1,m2,m3,m4,m5,m6,m7,m8,m9,\
        m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,\
        m26,m1_err,m2_err,m3_err,m4_err,m5_err,m6_err,m7_err,m8_err,\
        m9_err,m10_err,m11_err,m12_err,m13_err,m14_err,m15_err,m16_err,\
        m17_err,m18_err,m19_err,m20_err,m21_err,m22_err,m23_err,m24_err,\
        m25_err,m26_err = np.loadtxt(phot_name[ff]+'.trimmed',unpack=True)

        # Record various properties of observation (chip,amp,mjd, etc.).
        if ext == 0:   # LACOSMIC cleaned
            dd=pyfits.open(file_name)
            hdr0=dd[0].header  
            im = dd[0].data    
            dd.close()
            amps = hdr0['CCDAMP']
            filt = hdr0['FILTER']
            LTV1 = hdr0['LTV1']
            LTV2 = hdr0['LTV2']
            chip = hdr0['CCDCHIP']  # Change for driz data, 2
            shut = hdr0['SHUTRPOS']  # Change for driz data, eg, 'A'
            expt = hdr0['EXPTIME']
            mjd_avg = (hdr0['EXPEND'] - hdr0['EXPSTART'])/2. + hdr0['EXPSTART']
            # time between observation starts in minutes
            mjd_deltat = (hdr0['EXPEND'] - hdr0['EXPSTART'])*24.0*60.0                
            xcp = xc - LTV1                #physical pixel location
            ycp = yc - LTV2
            sizaxis1 = hdr0['SIZAXIS1']
            sizaxis2 = hdr0['SIZAXIS2']
            
        elif ext == 1:
            dd=pyfits.open(file_name)
            hdr0=dd[0].header 
            hdr1=dd[1].header 
            im = dd[1].data    
            dd.close()
            amps = hdr0['CCDAMP']
            filt = hdr0['FILTER']
            LTV1 = hdr1['LTV1']
            LTV2 = hdr1['LTV2']
            chip = hdr1['CCDCHIP']  
            shut = hdr0['SHUTRPOS'] 
            expt = hdr0['EXPTIME']
            mjd_avg = (hdr0['EXPEND'] - hdr0['EXPSTART'])/2. + hdr0['EXPSTART']
            # time between observation starts in minutes
            mjd_deltat = (hdr0['EXPEND'] - hdr0['EXPSTART'])*24.0*60.0                
            xcp = xc - LTV1                #physical pixel location
            ycp = yc - LTV2
            sizaxis1 = hdr1['SIZAXIS1']
            sizaxis2 = hdr1['SIZAXIS2']
                               
        print "...UPDATING DICTIONARY..."
        data[ff] = {'filename':file_name, 'amp':amps,'shutter':shut,\
                    'mjd_avg':mjd_avg, 'mjd_deltat': mjd_deltat, \
                    'chip': chip, 'axis1':sizaxis1, 'axis2':sizaxis2, \
                    'xc':xc, 'yc':yc,'xcp':xcp, 'ycp':ycp, \
                    'background': back, 'background_rms':backrms, \
                    'exptime': expt, \
                    'flux':[f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,\
                            f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,\
                            f25,f26], \
                    'mag':[m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,\
                            m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,\
                            m25,m26], \
                    'merr':[m1_err,m2_err,m3_err,m4_err,m5_err,m6_err,\
                            m7_err,m8_err,m9_err,m10_err,m11_err,m12_err,\
                            m13_err,m14_err,m15_err,m16_err,m17_err,\
                            m18_err,m19_err,m20_err,m21_err,m22_err,\
                            m23_err,m24_err,m25_err,m26_err]}
            
        """if images_seperate:
            # Construct diagnostic diagrams.
            make_source_location_histogram_plots_uvis(data, file_name, ff, im, \
                                                        find_name[ff], \
                                                        filt, \
                                                        path_to_photcat+\
                                                        'temp_lacos/')
            
            # Remove all unneeded temporary files.
            tmp = glob.glob(path_to_photcat+'temp_lacos/*cnts.fits.back*') + \
                  glob.glob(path_to_photcat+'temp_lacos/*cnts.fits') + \
                  glob.glob(path_to_photcat+'temp_lacos/*.trimmed*')
            for tt in tmp: os.remove(tt)
        else:
            # Construct diagnostic diagrams.
            make_source_location_histogram_plots_uvis(data, file_name, ff, \
                                                        im, find_name[ff], \
                                                        filt, \
                                                        path_to_photcat)
            
            # Remove all unneeded temporary files.
            tmp = glob.glob(path_to_photcat+'*cnts.fits.back*') + \
                  glob.glob(path_to_photcat+'*cnts.fits') + \
                  glob.glob(path_to_photcat+'*.trimmed*')
            for tt in tmp: os.remove(tt)"""
    
        # Save ascii catalog of sources.
        print "...WRITING CATALOG..."
        # Note we write to catalog for every iteration so that we don't 
        # loose work if we are in manual mode and DS9 decides to throw a tantrum.
        make_photom_catalog_uvis(data, filt, path_to_photcat)
    
    # Append the new photcat file to the old photcat file (if it exists). 
    # Then delete the new file and rename the old back
    # to <filt>_photcat.dat.
    append_to_photcat=True
    if append_to_photcat==True:
		photcatstore_list = glob.glob(path_to_photcat + '*photcat.dat')
		if photcatstore_list != []:
			append_phot_catalog(filt, perm_path_to_photcat,path_to_photcat)
			photcat_list = glob.glob(perm_path_to_photcat + '*photcat.dat')
			for photcat in photcatstore_list:
				os.remove(photcat)
        #rename_photcat(filt, origin=perm_path_to_photcat, revert=True)
    if append_to_photcat==False:
    	print '!!!!!!!!!!!!!!!!!!!!!!!!!!'
    	print 'not appending to photcat. data will be in seperate catalog'
    
    # Move all files out of temporary `temp_lacos`, if needed.
    
    
    # Move the position-histogram plots.
    #move_files(origin=path_to_photcat)    
    print "PTSRC_PHOTOM_UVIS COMPLETE"
    
    
#-------------------------------------------------------------------------------#
    
def main_photom():
	#set these first!
	data_dir=set_paths()['data']
	objects = set_paths()['objects']
	temp_for_lacosmicc=True
	append_to_photcatt=True #if true, will append new phot to existing catalog. otherwise, 
						   #the data will stay in its own <filter>_photcat.dat in temp directory
	
	move_flts=True #set this True if you wish to move files out of temp_for_lacosmic folder
				   #to the main data storage directory after being processed.
	
	for object in objects:
		filters=set_paths()['filters']

		if temp_for_lacosmicc:
			paths=[data_dir+'/' + object + '/'+ f + '/temp_for_lacosmic/flt_cleans/' for f in filters]
		else: #not grabbing flts from temp directory, using ones in permanent flt_cleans location
			paths=[data_dir+'/' + object + '/'+ f + '/flt_cleans/' for f in filters]
			
		perm_paths_to_photocat=[data_dir+'/'+object+'/'+f+'/flt_cleans' for f in filters]
		path_to_fltss=[]
		perm_path_to_photcatss=[]
		
		for i, path_to_flts in enumerate(paths):
			path_to_fltss.append(path_to_flts)
			perm_path_to_photcatss.append(perm_paths_to_photocat[i])
			
		
			if not os.path.isdir(perm_paths_to_photocat[i]):
				print 'no existing photcat found, making directory'
				os.makedirs(perm_paths_to_photocat[i])
			
			
			if len(filters) <= 8:
				num_cores=len(filters)
			else:
				num_cores=8
		
		args=zip(path_to_fltss,perm_path_to_photcatss,[append_to_photcatt]*len(path_to_fltss))
		
				
		Pool(3).map(ptsrc_photom_uvis, args)
		
		if temp_for_lacosmicc and move_flts:
			filters=set_paths()['filters']
			move_files_from_temp_directories([object],filters)

if __name__=='__main__':
	main_photom()

