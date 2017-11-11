"""Runs ``calwf3`` in current working directory.
    
Author
------

    C. Shanahan
    
Use
---

   cd into directory containing RAWS. Then,
    
    >>> python run_calwf3.py
    
    Can also be imported by other functions and supplied a path to location of files.
    
Outputs
-------

    FLT files. Removes RAW files.
    

"""

import glob
import os
import shutil
from wfc3tools import calwf3
from astropy.io import fits

def change_pctecorr_keyword(file_path):
    """Changes the 'PCTECORR' keyword in the RAW file to 'OMIT'.
    
    Parameters:
        file_path : string
            Path to the RAW FITS file.
            
    Returns:
        nothing
        
    Outputs:
        original RAW file with the changed header key.         
    """
    raw = fits.open(file_path, mode='update')
    hdr = raw[0].header
    
    if 'PCTECORR' in hdr.keys():
        hdr['PCTECORR'] = 'OMIT'
    
    raw.close()


def run_calwf3(path_to_files = ''):

    """Recalibrate with UNITCORR set to 'OMIT'.
    """
    raws = glob.glob(path_to_files+"/*raw.fits")
    print raws
    # Change the keywords.
    for raw in raws:
        change_pctecorr_keyword(raw)
    
    # Recalibrate.
    for raw in raws:
        calwf3(raw)
        
    # Remove TRA files and linearity corr file.    
    tras = glob.glob(path_to_files+"/*tra.fits")
    for tra in tras:
        os.remove(tra)

if __name__ == '__main__':
	run_calwf3()