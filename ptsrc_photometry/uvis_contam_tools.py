#! /usr/bin/env python

"""Contains often-used functions for UVIS contam-stability monitor.
    
Author:

    C.M. Gosmeyer, Aug. 2014
    
    D.M. Hammer, 2012

Use:

    From here import functions.

Notes:

    Rule of thumb: If I import the same function in three seperate
    and very different scripts, it is time to place it here.

"""
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import glob
import os
import shutil
import sys
import numpy as np
from astropy.io import fits
from contextlib import contextmanager

from photutils import CircularAperture, DAOStarFinder, aperture_photometry
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd 
import warnings

def fxn():
    warnings.warn("deprecated", np.VisibleDeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
    
def sub_amp():
    dictt={'UVIS1-C512A-SUB':(1,'A'),'UVIS1-C512B-SUB':(1,'B'),\
           'UVIS2-C512C-SUB':(2,'C'), 'UVIS2-C512D-SUB':(2,'D'),\
           'UVIS2-M512C-SUB':(2,'C')}
    return dictt

def make_PAMcorr_image(image, outfile='default', extension=0):
    """Creates the Pixel Area Map (PAM) image.
    
    Parameters:
        image : string
            Name of FITS file.
        outfile : string
            Name of outfile.
        extension : int 
            Extension=0 (UV images) by default.
            
    Returns:
        outfile : FITS file
            The pixel area map called ``<image rootname>_PAM.fits``
            
    Outputs:
        nothing
    """
    
    pamdir = '/grp/hst/wfc3t/cshanahan/'
    
    # -- Parse output filename & save a copy to file (NOTE: if outfile == input image, data is overwritten).
    if (image != outfile):
        if outfile == 'default': 
            outfile = image.split('.fits')[0] + '_PAM.fits'
        shutil.copy(image,outfile)
    else:
        print 'OVERWRITING DATA FOR IMAGE: '+image+'.'
    
    # -- Read in fits image and header (assume flt/flc - should handle both full- & sub-arrays)
    prihdr = fits.getheader(outfile)
    fdata = fits.getdata(outfile, ext=extension)   # ext=1 for IR?
    exptime = prihdr['exptime']
    detector = prihdr['detector']
    
    # -- Cycle through each SCI extension
    hdulist = fits.open(outfile,mode='update')
    for ff in xrange(len(hdulist)):
        if hdulist[ff].name == 'SCI':
            
            # -- read in header and data info
            scihdr = hdulist[ff].header
            data = hdulist[ff].data
            if detector == 'IR':
                chip = 1
            elif detector == 'UVIS':
                chip = scihdr['CCDCHIP']
            else:
                raise Exception('Detector '+detector+' not covered in our case list.')
            
            naxis1 = scihdr['NAXIS1']
            naxis2 = scihdr['NAXIS2']
            x0 = int(np.abs(scihdr['LTV1']))
            y0 = int(np.abs(scihdr['LTV2']))
            x1 = int(x0 + naxis1)
            y1 = int(y0 + naxis2)
            
            # -- apply the PAM
            if detector == 'UVIS':
                if chip == 1:
                    pam=fits.getdata(pamdir+'UVIS1wfc3_map.fits')
                    
                    hdulist[ff].data = data * pam[y0:y1,x0:x1]

                elif chip == 2:
                    pam=fits.getdata(pamdir+'UVIS2wfc3_map.fits')
                    
                    hdulist[ff].data = data * pam[y0:y1,x0:x1]

                else:
                    raise Exception('Chip case not handled.')
                    hdulist.close()
            elif detector == 'IR':
                pam=fits.getdata(pamdir+'ir_wfc3_map.fits')
                hdulist[ff].data = data * pam[y0:y1,x0:x1]
            else:
                raise Exception('Detector '+detector+' not covered in our case list.')
                hdulist.close()
    hdulist.close()

    return outfile


#-------------------------------------------------------------------------------#

def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=0, return_array=0, return_median=0):
    """Computes an iteratively sigma-clipped mean on a data set. 
    
    Clipping is done about median, but mean is returned by default 
    (use return_median to return the clipped median value).
    
    Parameters:
        indata : array_like
            Input data.
        
        clipsig : float
            Number of sigma at which to clip.
        
        maxiter : int
            Ceiling on number of clipping iterations.
        
        converge_num : float
            If the proportion of rejected pixels is less than
            this fraction, the iterations stop.
        
        verbose : {0, 1}
            Print messages to screen?
        
        return_array : {0, 1}
            Return the final array indices that were used to compute
            statistics.
            
        return_median : {0, 1}
            Return the median if 1. Return the mean if 0. 
        
    Returns:
        val : float
            The N-sigma clipped mean or median.
        sigma : float
            The standard deviation of remaining pixels.
        
    Outputs:
        Prints to screen mean or median statistics.      
        
    History:
        * 21/10/1998 Written by RSH, RITSS
        * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
        * 24/11/2009 Converted to Python. PLL.
        * 08/01/2013 Added option to return the array indices of non-clipped pixels. DMH
        * Added option to return median of the clipped array. DMH
    
    Notes:
        This is based on MYMEANCLIP routine from ACS library.    
    
    Examples:
    
    >>> mean, sigma = meanclip(indata)
    """
    # Flatten array
    skpix = indata.reshape( indata.size, )
    
    # initialize array to store indices of pixels used to compute stats
    arrind = np.arange(0,skpix.size)
    
    ct = indata.size
    iter = 0; c1 = 1.0 ; c2 = 0.0
    
    while (c1 >= c2) and (iter < maxiter):
        lastct = ct
        medval = np.median(skpix)
        sig = np.std(skpix)
        wsm = np.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]
            arrind = arrind[wsm]
        
        c1 = abs(ct - lastct)
        c2 = converge_num * lastct
        iter += 1
    
    if return_median:
        val = np.median(skpix)
        val_type = 'median'
    else:
        val  = np.mean( skpix )
        val_type = 'mean'
    sigma = np.std( skpix )
    
    if verbose:
        if return_median:
            prf = 'MEDIANCLIP:'
            print '%s %.1f-sigma clipped median' % (prf, clipsig)
            print '%s Median computed in %i iterations' % (prf, iter)
            print '%s Median = %.6f, sigma = %.6f' % (prf, val, sigma)
        
        else:
            prf = 'MEANCLIP:'
            print '%s %.1f-sigma clipped mean' % (prf, clipsig)
            print '%s Mean computed in %i iterations' % (prf, iter)
            print '%s Mean = %.6f, sigma = %.6f' % (prf, val, sigma)
    
    if return_array:
        return np.copy(arrind)
    else:
        return val, sigma

def circular_mask(arr_shape, r, x_offset=0, y_offset=0):
    """Generates circular mask for 2D image. Function from the IRAF photometry script. 
        
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
    
def diagnostic_source_finding_plots(ifile,coo_tab=None):

    hdu = fits.open(ifile)
    data = hdu[0].data

    norm = ImageNormalize(stretch=LogStretch())
    plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
    
    if coo_tab:
        positions = (np.array(coo_tab['xcentroid'].tolist()), \
                    np.array(coo_tab['ycentroid'].tolist()))
        apertures = CircularAperture(positions, r=10.)
        apertures.plot(color='blue', lw=2, alpha=1)
    
    plt.show()
    
def get_wfc3_zeropoint(filter):
    """Gets (infinite aperture) zeropoint for WFC3 filters. 
    
    Parameters:
        filter : string
            Name of the filter.
    
    Returns:
        zp[filter.upper()] : float
            Zeropoint of the filter.
    
    Outputs:
        nothing
    
    Notes:
        Pre two-chip solution.
        Add filters as needed.
        
    References:
        Infinite aperture ABmag column http://www.stsci.edu/hst/wfc3/phot_zp_lbn
    """
    # array of WFC3 AB zeropoints (not all here - add as needed)
    zp = {'F098M':25.6674, 'F105W':26.2687, 'F110W':26.8223, 'F125W':26.2303, \
          'F126N':22.8609, 'F127M':24.6412, 'F128N':22.9726, 'F130N':22.9900, \
          'F132N':22.9472, 'F139M':24.4793, 'F140W':26.4524, 'F153M':24.4635, \
          'F160W':25.9463, 'F164N':22.9089, 'F167N':22.9568, 'F200LP':27.3407, \
          'F218W':22.9641, 'F225W':24.0403, 'F275W':24.1305, 'F280N':20.9536, \
          'F300X':24.9565, 'F336W':24.6682, 'F343N':23.8916, 'F350LP':26.9435, \
          'F373N':21.9052, 'F390M':23.6360, 'F390W':25.3562, 'F395N':22.6646, \
          'F410M':23.5962, 'F438W':24.8206, 'F475X':26.1579, 'F467M':23.6799, \
          'F469N':21.7911, 'F475W':25.6799, 'F487N':22.2174, 'F502N':22.2956, \
          'F547M':24.7467, 'F555W':25.7906, 'F600LP':25.8746, 'F606W':26.0691, \
          'F621M':24.6001, 'F625W':25.5263, 'F631N':21.8922, 'F645N':22.2251, \
          'F656N':20.4827, 'F657N':22.6408, 'F658N':21.0287, 'F665N':22.4899, \
          'F673N':22.5684, 'F680N':23.7992, 'F689M':24.4641, 'F763M':24.2070, \
          'F775W':24.8555, 'F814W':25.0985, 'F845M':23.7811, 'F850LP':23.8338, \
          'F953N':20.3720}
    
    
    if zp.has_key(filter.upper()): 
        return zp[filter.upper()]
    else: 
        raise Exception('Zeropoint is not specified for this filter: '+filter)
        
def get_UVIS_gain(quad):

    
    gain = {'A':1.56,'B':1.55,'C':1.58,'D':1.57}
    
    return gain[quad]
    
