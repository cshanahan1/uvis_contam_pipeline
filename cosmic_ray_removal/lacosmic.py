#! /usr/bin/env python
"""
	This script will run the Python ro

"""

import shutil
import glob
import os
import sys
sys.path.append("/grp/hst/wfc3t/cshanahan/uvis_contam_monitor/uvis_contam_scripts")
from set_paths_and_params import *
from lacosmic_param_dictionary import lacosmic_param_dictionary
import pylab
from astropy.io import fits as pyfits
import cosmics
import img_scale

def create_diagnostic_png(filename, outfilename='Default'):

    pylab.ioff()
    # create page for plots
    page_width = 21.59/2
    page_height = 27.94/2
    fig = pylab.figure(figsize=(page_width, page_height))

    file_clean = (filename.split('.fits')[0]+'.clean.fits')
    file_mask = (filename.split('.fits')[0]+'.mask.fits')

    scmax = 7000 # scale_max for raw and clean
    scmin = 3 # scale_min for raw and clean

    # Plot the original image
    pylab.subplot(3,2,1, aspect='equal') # 311
    image_orig = pyfits.open(filename)
    image_orig_ext = image_orig[1].data
    image_orig_scaled = img_scale.log(image_orig_ext, \
                                      scale_min=scmin, \
                                      scale_max=scmax)
    plt_orig = pylab.imshow(image_orig_scaled, aspect='equal')
    pylab.title('Original (SCI)')

    # Plot cut of original image
    pylab.subplot(3,2,2, aspect='equal')
    image_orig_cut = img_scale.log(image_orig_ext[175:275,175:275], \
                                   scale_min=scmin, \
                                   scale_max=scmax)
    plt_orig_cut = pylab.imshow(image_orig_cut, aspect='equal')
    pylab.title('Original (SCI)')
    image_orig.close()


    # Plot the mask image
    pylab.subplot(3,2,3) #312
    image_mask = pyfits.open(file_mask)
    image_mask_ext = image_mask[0].data
    plt_orig = pylab.imshow(image_mask_ext, aspect='equal', vmin=-2, vmax=1) #-2, -5
    pylab.title('Mask')

    # Plot cut of the mask image
    pylab.subplot(3,2,4)
    plt_mask_cut = pylab.imshow(image_mask_ext[175:275,175:275], \
                                aspect='equal', vmin=-2, vmax=1)
    pylab.title('Mask')
    image_mask.close()


    # Plot the LACosmic-cleaned image
    pylab.subplot(3,2,5) #313
    image_clean = pyfits.open(file_clean)
    image_clean_ext = image_clean[0].data
    image_clean_scaled = img_scale.log(image_clean_ext, \
                                       scale_min=scmin, \
                                       scale_max=scmax)
    plt_clean = pylab.imshow(image_clean_scaled, aspect='equal')
    pylab.title('Clean')

    #Plot cut of the LACosmic-cleaned image
    pylab.subplot(3,2,6)
    image_clean_cut = img_scale.log(image_clean_ext[175:275,175:275], \
                                    scale_min=scmin, scale_max=scmax)
    plt_clean_cut = pylab.imshow(image_clean_cut, aspect='equal')
    pylab.title('Clean')
    image_clean.close()

    if outfilename == 'Default':
        pylab.savefig(filename.split('.fits')[0] + '.png')
    else:
        pylab.savefig(outfilename)
    pylab.close()
    pylab.ion()


	
def sort_files(origin,temp_folder,create_pngs=False):
    """Moves clean and mask FITS images and their PNGs into their own
    subdirectories, `flt_cleans`., `flt_masks`, and `png_masks_cleans`.

    Parameters:
        origin : string
            Path to where the mask and clean FLTs and the PNGs are located.
            If left blank, assume Current Working Directory.
        dest : string
            Path where you want the directories to go. (Usually will want
            origin and dest to match). If left blank, assume
            Current Working Directory.
        temp_folder : {True, False}
            False by default. Switch on if want the clean files placed
            in `flt_cleans/temp/`.


    Returns:
        nothing

    Outputs:
        nothing
    """
    
    cleans=glob.glob(origin+'*.clean.fits')
    masks=glob.glob(origin+'*.mask.fits')
   
    
   
    folder_cleans=origin+'flt_cleans'
    folder_masks=origin+'flt_masks'
    	
    	
    if not os.path.isdir(folder_cleans):
    	os.makedirs(folder_cleans)
    for clean in cleans:
    	if not os.path.isfile(folder_cleans+'/'+clean.replace(origin,'')):
    		shutil.move(clean,folder_cleans)
    	else:
    		print clean.replace(origin,'') + ' already exists in flt_cleans'
    		
    if not os.path.isdir(folder_masks):
    	os.makedirs(folder_masks)
    	
    for mask in masks:
    	if not os.path.isfile(folder_masks+'/'+mask.replace(origin,'')):
    		shutil.move(mask,folder_masks)
    	else:
    		print mask.replace(origin,'') + ' already exists in flt_masks'
    		
    if create_pngs:
    	folder_pngs=origin+'/'+'pngs'
    	pngs=glob.glob(origin+'/'+'*.png')
    	if not os.path.isdir(folder_pngs):
    		os.makedirs(folder_pngs)
    	for png in pngs:
    		shutil.move(png,folder_pngs) 
    			
    		
    


def run_lacosmic(filename, sigclip, sigfrac, objlim, niter, sigclip_pf):
	
	filename=str(filename)
	sigclip = float(sigclip)
	sigfrac = float(sigfrac)
	objlim = int(objlim)
	niter = int(niter)
	sigclip_pf = float(sigclip_pf)
    

	fits_file = pyfits.open(filename)
	arrayy=fits_file[1].data
	headerr=fits_file[1].header+fits_file[0].header
    
	flshcorr = headerr['FLSHCORR']
	fits_file.close()

	if sigclip_pf == 0.0 or flshcorr == 'OMIT':
		sigclip = sigclip
		print 'FLSHCORR set to OMIT.'
	elif flshcorr == 'COMPLETE':
		sigclip = sigclip_pf
		print 'FLSHCORR set to COMPLETE.'
    
	# Build the object :
	c = cosmics.cosmicsimage(arrayy, gain=1.5, readnoise=3.0, sigclip = sigclip, sigfrac = sigfrac, objlim = objlim, verbose=False)
	# cosmics module is very verbose, if you wish to turn off check out cosmics.py
	c.run(maxiter = niter, verbose=False)
	print 'DONE'
	
	cosmics.tofits(filename.replace('flt.fits','flt.clean.fits'), c.cleanarray, headerr, verbose=False)
	
	cosmics.tofits(filename.replace('flt.fits','flt.mask.fits'), c.mask, headerr,verbose=False)


def run_lacosmic_main(origin,temp_folder, create_pngs):
	print "ORIGIN IS " + str(origin)
	param_dict = lacosmic_param_dictionary()
	#fits_list = glob.glob(origin + '/'+ '*fl*.fits')
	fits_list=glob.glob(origin + '/'+ '*flt.fits')
	
	filt = pyfits.open(fits_list[0])[0].header['FILTER']
  

    # Run LACOSMIC.
	for fitss in fits_list:
		if filt not in param_dict.keys():
			print "Filter not in Param Dictionary. Using default values."
			if 'N' in filt:
				sigclip = 4.5
				objlim = 5
			else:
				sigclip = 5.0
				objlim = 2
			sigfrac = 0.3
			niter = 3
			sigclip_pf = 9.5

		else:
			sigclip = param_dict[filt][0]
			sigfrac = param_dict[filt][1]
			objlim = param_dict[filt][2]
			niter = param_dict[filt][3]
			sigclip_pf = param_dict[filt][4]

		run_lacosmic(fitss, sigclip, sigfrac, objlim, niter, sigclip_pf)
		if create_pngs:
			create_diagnostic_png(fitss)
	
	sort_files(origin=origin, temp_folder=temp_folder, create_pngs=create_pngs)


if __name__ == '__main__':

    try:
    	flt_dir = str(sys.argv[1])
    except:
    	print "usage: run_lacosmic.py path/to/directories/with/flts"
    
	"""by default, output will be in directory with flts.
	if temp_folder, a flt_cleans/temp will be created for photometry.
	otherwise, the output will be in a subdiectory in the flts directory
	called flt_cleans"""
    run_lacosmic_main(origin=flt_dir, temp_folder=True, create_pngs=False)

    print "Finished at last."
	
	
	
	
	