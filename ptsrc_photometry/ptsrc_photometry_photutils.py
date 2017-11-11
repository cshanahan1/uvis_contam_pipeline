from astropy.io import ascii, fits
from astropy.stats import mad_std
import copy
import glob
from photutils import CircularAperture, DAOStarFinder, aperture_photometry
import numpy as np
import os.path
import pandas as pd 
import shutil
from uvis_contam_tools import circular_mask, get_UVIS_gain, get_wfc3_zeropoint, meanclip,\
     make_PAMcorr_image, diagnostic_source_finding_plots, sub_amp
from multiprocessing import Pool
import warnings
from set_paths_and_params import *
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*    
def make_pamcorr(ifile):
    
    """ Makes PAM corrected image"""
    
    make_PAMcorr_image(ifile,outfile=ifile+'.pamcorr', extension=0)
    ifile = ifile+'.pamcorr'
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*       
def create_photometry_file(phot_file_path,aperture_sizes):

    """ Creates the photometry file and writes out the header with column names, if 
        the file does not exist yet. 
    
     Parameters
     ----------
        phot_file_path : str
            Full path to output photometry file.
        aperture_sizes : list of ints
            List of all aperture sizes that are used in the photometry.
     Outputs
     ------- 
        Creates the file at phot_file_path and writes the header to the first line.
     
     """


    flux_colnames = ['f{0},f{0}err'.format(str(ap_size)) for ap_size in aperture_sizes]
    flux_colnames = ','.join(flux_colnames)
    
    mag_colnames = ['m{0},m{0}err'.format(str(ap_size)) for ap_size in aperture_sizes]
    mag_colnames = ','.join(mag_colnames)
    
    phot_header = '#filename,prop_id,amp,shutter,mjd,chip,axis1,axis2,exptime,xc,yc,back,backrms,'
    phot_header = phot_header + flux_colnames + ','+mag_colnames + '\n'
    
    print phot_header
    with open(phot_file_path,'w') as f:
        f.write(phot_header)
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*        
def fetch_header_data(ifile,hdu,paths):

    """ Gets parameters from the header of the file to write out to the photometry table.
    
     Parameters
     ----------
        ifile : str
            Path of input file.
        hdu : hdu list
            HDU list of input file.
     Returns
     ------- 
        
        header_info : list
            List of info from header, [filename,prop_id,amp,shutter,mjd,chip,axis1,axis2,exptime].
     
     
     """
     
    filename = os.path.basename(ifile)[0:9]
    header = hdu[0].header 
    
    subarray_amp=sub_amp()
    
 
    subarray = header['APERTURE']
    
    prop_id = header['PROPOSID']

    #print 'subarray ' , subarray
    try:
        chip, amp = subarray_amp[subarray]
    except KeyError:
        print 'observaion not in desired subarray. moving to bad_files/wrong_subarray' 
        bad_files_dir = paths['data']+'/bad_files'
        if not os.path.isdir(bad_files_dir):
            os.mkdir(bad_files_dir)
        if not os.path.isdir(bad_files_dir+'/wrong_subarray'):
            os.mkdir(bad_files_dir+'/wrong_subarray')
        #also move the flt.fits and mask.fits
        data_dir = paths['data']+'/data'
        flts = glob.glob(data_dir+'/{}*'.format(filename))
        flts = flts + glob.glob(data_dir+'temp_for_photometry/{}*.fits'.format(filename))
        masks = glob.glob(data_dir+'/flt_masks/{}*.fits'.format(filename))
        masks=masks+glob.glob(data_dir+'temp_for_photometry/flt_masks/{}*.fits'.format(filename))
        move_files = flts + masks + glob.glob(ifile.replace(os.path.basename(ifile),filename+'*.fits'))

        for f in move_files:
            shutil.move(f, paths['data']+'/bad_files/wrong_subarray')
        #os.remove(ifile)
        print 'removed '  + ifile
        
        return 0 
     
    shutter = header['SHUTRPOS']
     
    mjd = header['EXPSTART']
    
    exptime = header['EXPTIME']
    
    axis1=header['NAXIS1']
    axis2=header['NAXIS2']
    
    header_info = [filename,prop_id,amp,shutter,mjd,chip,axis1,axis2,exptime]
    
    

    return header_info
    
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*     
def run_dao_star_finder(ifile,image_name,hdu, write_coo_tab = True, 
                        bkgrnd_threshold = True, fwhm = 2.5, threshold = 1000):
                        
    """ Runs photutils.DAOStarFinder on input fits file and finds the single source. 
        If multiple sources / no sources are found on an observation, they are moved to a 
        sub directory called 'failed_source_finding', that is created in the same 
        directory as the input file if it doesn't exist. Returns an astropy table with
        source x and y position. 
    
     Parameters
     ----------
        ifile : int
            Input file path.
        exclude_border : bool
            True if you wish to exclude a 20 pixel border.

     Returns
     -------
     
        coo_tab : astropy table 
            Table with x and y positions of source
            
        If no source/multiple sources are found, this function will return False and 
        move the offending file to the 'failed_source_finding', a subdirectory of the 
        input file location. 
        
     """
     
    data = hdu[0].data
    
    fname = os.path.basename(ifile)[0:9]

    if bkgrnd_threshold:
        threshold = threshold * mad_std(data)

    daofindd = DAOStarFinder(threshold=threshold, fwhm=fwhm, exclude_border = True)
    coo_tab = daofindd(data)
    
    if len(coo_tab) != 1:
        if not os.path.isdir(ifile.replace(os.path.basename(ifile),'failed_source_finding')):
            os.mkdir(ifile.replace(os.path.basename(ifile),'failed_source_finding'))
            os.mkdir(ifile.replace(os.path.basename(ifile),'failed_source_finding/none'))
            os.mkdir(ifile.replace(os.path.basename(ifile),'failed_source_finding/multiple'))

    if len(coo_tab) < 1:
        print 'No sources found for'+ os.path.basename(ifile) +'. Moving to next file.' 
        
        move_files = glob.glob(ifile.replace(os.path.basename(ifile),'')+fname+'*.fits')
        
        #also move the flt, one directory back from the current flt_clean ifile
        
        flt_dir = '/'.join(ifile.split('/')[:-2])
        move_files = move_files + glob.glob(flt_dir+'/{}*'.format(fname))
        
        print 'move_files', move_files
        for filee in move_files:
            shutil.move(filee, ifile.replace(os.path.basename(ifile),'failed_source_finding/none'))
            
        #os.remove(ifile)
        return False
        
    if len(coo_tab) > 1:
    #else:
        print 'Multiple sources found for'+ os.path.basename(ifile) +'. Moving to next file.' 
    
        move_files = glob.glob(ifile.replace(os.path.basename(ifile),'')+fname+'*.fits')
        
        #also move the flt, one directory back from the current flt_clean ifile
        
        flt_dir = '/'.join(ifile.split('/')[:-2])
        move_files = move_files + glob.glob(flt_dir+'/{}*'.format(fname))

        print 'move_files', move_files
        for filee in move_files:
            shutil.move(filee, ifile.replace(os.path.basename(ifile),'failed_source_finding/multiple'))
        #os.remove(ifile)
        return False 
         
    return coo_tab
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*                  
def calculate_bkgrnd(hdu, filename, coo_tab, mask_rad = 80,mask_border_size=20):

    """ Calculates background and background RMS for an image. Uses the meanclip function 
        for an iterative sigma clipping of the sky region. Sky level is determined from 
        image pixels excluding circular region around the source, and a border excluded 
        around the edges.
        
        Parameters
        ----------
            hdu: HDU list
                HDU list from input file
            filename : string
                Basename of input file.
            coo_tab : astropy table

        Returns
        -------
            back : float
                Background sky level returned by the iterative sigma clipping routine 
                meanclip on non-excluded regions.
            backrms : float
                Background RMS.
    """
   
   
    data = hdu[0].data 
    
    xc = coo_tab['xcentroid'][0]
    yc = coo_tab['ycentroid'][0]
    
    naxis1=hdu[0].header['NAXIS1']
    naxis2=hdu[0].header['NAXIS2']
    
    #create a temporary image to mask excluded zones
    temp_data = copy.copy(data)
    
    #mask region around source
  
    temp_data[circular_mask(temp_data.shape, mask_rad, x_offset=(xc-naxis1/2.0), 
    y_offset=(yc-naxis2/2.0))] = -99999.0

    #mask pixels around border
    temp_data[:,0:mask_border_size] = -99999.0
    temp_data[:,-mask_border_size:] = -99999.0
    temp_data[0:mask_border_size,:] = -99999.0
    temp_data[-mask_border_size:,:] = -99999.0
    
    #now run meanclip.py on the masked data. 
    
    #generate initial guess for lower/upper limits (use 10 sigma).
    fmasked_data = np.ndarray.flatten(temp_data)
    llim = -100
    ulim = 10000.0
    init_median,init_rms = meanclip(fmasked_data[(fmasked_data > llim) & \
                                    (fmasked_data < ulim)], maxiter=7, \
                                    return_median=1)
    llim = init_median - 10.0*init_rms
    ulim = init_median + 10.0*init_rms
    
    #calculate backgound and rms
    back,backrms=meanclip(fmasked_data[(fmasked_data > llim) & (fmasked_data < ulim)], \
                          maxiter=7)
                          
    return (back, backrms)
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*  
def calculate_mag_errs(flux, mag, aperture_size, back, backrms, epdau):
    """Calculate measurement errors in the same manner as IRAF/DAOphot.
    
    Parameters
    ----------
    flux : float
        Sky subtracted flux of star.
    back : float
        Sky level.  
    backrms : float
        Background RMS.
    epdau : float
        CCD gain
     
     Returns
     --------
     errs : tuple
        (error in insturmental magnitudes, error in flux)
        
    """
    
    #convert e- to ADU
    flux = flux / epdau 
    back = back / epdau 
    backrms = backrms / epdau
    
    
    area = np.pi * (aperture_size/2)**2
    
    err_terms = (flux/epdau) + (area*backrms**2) + (area*backrms**2/back)
    
    flux_err = np.sqrt(np.abs(err_terms)) # in ADU
    
    flux_err = flux_err * epdau
    
    merr = 1.0857*flux_err/(flux * epdau)
    
    return (merr,flux_err)
            
    
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*   
def run_aperture_photometry(ifile, hdu, coo_tab, back, backrms, aperture,PAMcorr=True):

    """Runs photutils.aperture_photometry on an observation, for a single circular
        aperture size. Subtracts the (given) sky level from the aperture."""

    header = hdu[0].header
    exptime = header['EXPTIME']
    filterr=header['FILTER']
    data = hdu[0].data 
    
    
    xc = coo_tab['xcentroid'][0]
    yc = coo_tab['ycentroid'][0]

    apertures = CircularAperture((coo_tab['xcentroid'], coo_tab['ycentroid']), r=aperture)
    phot_table = aperture_photometry(data - back, apertures)

        
    flux = phot_table['aperture_sum'][0]
    mag = -2.5*np.log10(flux/exptime)+get_wfc3_zeropoint(filterr)
    
    subarray_amp=sub_amp()
     
    subarray = header['APERTURE']
    chip, amp = subarray_amp[subarray]
    
    epdau = get_UVIS_gain(amp)
    
    mag_err, flux_err = calculate_mag_errs(flux, mag, aperture, back, backrms, epdau)

    return (flux,flux_err,mag,mag_err)
    
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*  
def run_ptsrc_photometry(ifile, PAMcorr = True):

    #for i, ifile in enumerate(ifiles):
    
    paths = set_paths_and_params()
    
    try:
    	hdu_list = fits.open(ifile)
    except IOError:
    	print 'fits file corrupted. removing file'
    	os.remove(ifile)
    	return 0
    	
    
    header = hdu_list[0].header
    sci = hdu_list[0].data
    
    filename = os.path.basename(ifile)
    
    fits_data = (hdu_list, filename)
    
    filename = os.path.basename(ifile)[0:9]
    
    aperture_sizes = set_paths_and_params()['aperture_sizes'] 

    obs_info = fetch_header_data(ifile,hdu_list,paths)
    
    if not obs_info:
        return 0
        
    amp = obs_info[1]
 
    coo_tab = run_dao_star_finder(ifile,filename,hdu_list)
    
    if coo_tab: #if a single source was found
        xc,yc = coo_tab['xcentroid'][0],coo_tab['ycentroid'][0]
        obs_info.append(xc)
        obs_info.append(yc)
        back,backrms = calculate_bkgrnd(hdu_list, filename, coo_tab) 
        obs_info.append(back)
        obs_info.append(backrms)
        
        if PAMcorr:
            hdu_list.close()
            make_PAMcorr_image(ifile,outfile=ifile+'.pamcorr', extension=0)
            ifile = ifile+'.pamcorr'
            hdu_list = fits.open(ifile)
            header = hdu_list[0].header
            sci = hdu_list[0].data
        
        fluxes, flux_errs, mags, mag_errs = [],[],[],[]

        for ap_size in aperture_sizes:

            flux, flux_err,mag, merr = run_aperture_photometry(ifile, hdu_list, \
                                           coo_tab, back, backrms, ap_size)
    
            fluxes.append(flux)
            flux_errs.append(flux_err)
            
            mags.append(mag)
            mag_errs.append(merr)
            #print mag, merr
            
        for i, val in enumerate(fluxes):
            obs_info.append(fluxes[i])
            obs_info.append(flux_errs[i])
        for i, val in enumerate(mags):
            obs_info.append(mags[i])
            obs_info.append(mag_errs[i])
   
    else:
        hdu_list.close()    
        return 0
        
    hdu_list.close()
    #if PAMcorr:
        #os.remove(ifile)
        
    print 'done processing ' ,ifile
    
        
    return obs_info
    
    
    
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*
def main_ptsrc_photometry(ifiles, phot_file_path, paths, PAMcorr=True):
     
    aperture_sizes = paths['aperture_sizes']
              
    #create new photometry file if one doesnt exist
    if not os.path.isfile(phot_file_path):
        print 'file {} does not exist. creating'.format(os.path.basename(phot_file_path))
        create_photometry_file(phot_file_path,aperture_sizes)
                    

    #p = Pool(6)
    
    #Main returns information to be written to file rather than writing it in main.
    #Multiprocessing pools implement a queue to avoid collisions when subprocesses all
    #try to write to an output file, so the file should actually be written out here.

    
    """print 'Beginning photometry' 
    for result in p.map(run_ptsrc_photometry, ifiles):
        if result:
            obs_info_str = ','.join([str(item) for item in result])+'\n'
            with open(phot_file_path,'a') as f:
                f.write(obs_info_str)
                f.close()"""
                
    
  
    for ifile in ifiles:
        obs_info = run_ptsrc_photometry(ifile,paths)
        if obs_info:
            #print len(obs_info)
            with open(phot_file_path,'a') as f:
                obs_info_str = ','.join([str(item) for item in obs_info])+'\n'
                print obs_info_str
                f.write(obs_info_str)
                f.close()
    
    
    
#-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-*-*
if __name__ == '__main__':

    paths = set_paths_and_params()
    
    phot_file_path = 'F218W_photcat.dat'
    
    PAMcorr = True 
    
    ifiles = glob.glob('/grp/hst/wfc3t/cshanahan/UVIS_contam_monitor/data/data/GRW70/F218W/flt_cleans/*clean.fits')
    
    main_ptsrc_photometry(ifiles, phot_file_path, paths)
    


    
     